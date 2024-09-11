// Create your table
var table = new Tabulator("#resultsTable", {
    data: data, //assign data to table
    layout: "fitColumns", //fit columns to width of table
    columns: [ //Define Table Columns
        {title: "#", formatter:"rownum", width:40, headerSort:false},
        {title: "ID", field: "Unique_ID"},
        {title: "Name", field: "filename"},
        {title: "Type", field: "type"},
    ],
    pagination:"local",
    paginationSize:20,
    paginationSizeSelector:[10, 20, 50, 100],
});

table.on("cellClick", function(e, cell) {
    var value = cell.getValue();
    var field = cell.getField();

    if (field == "filename" || field == "Unique_ID") {
        rowData = cell.getRow().getData();
        document.getElementById('fileID').value = rowData["Unique_ID"];
        var enteredName = rowData["filename"];
        var enteredType = rowData["type"];

        create_query(rowData["Unique_ID"]);

        // Show the heading div after a delay so that data loading is complete
        setTimeout(function() {
            heading.classList.replace("d-none", "d-md-flex");
            searchField.classList.replace("d-md-flex", "d-none");
        }, 1000);

        var selectedName = document.getElementById('selectedName');
        selectedName.textContent = enteredName; 
        
    } else if (field == "Unit") {
        var url = "https://qudt.org/vocab/unit/" + encodeURIComponent(value);
        window.open(url, "_blank");
    } else if (field == "Value" && value == "Download") {
        var id = document.getElementById('fileID').value;
        fetch(`/rawdownload?id=${encodeURIComponent(id)}`, {
            method: 'GET'
        }).then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.blob();
        }).then(blob => {
            var url = window.URL.createObjectURL(blob);
            var a = document.createElement('a');
            a.href = url;
            a.download = `${id}.zip`;
            document.body.appendChild(a);
            a.click();
            a.remove();
        }).catch(error => {
            console.error('There has been a problem with your fetch operation:', error);
        });
    } else {
        console.log(`Geklickt auf Spalte: ${field} mit Value: ${value}`);
    }
});

// Creates a Query based on the specific input from the user and executes it
async function create_query(enteredName){
    
    // main logic for not extended queries
    try {
        const tableData = await executeSparqlQuery(enteredName);
        createTable(tableData);
    } catch (error) {
        document.getElementById('queryResults').textContent = error.message;
    }
}

// executes the sparql query and returns result as a promise
async function executeSparqlQuery(unique_id) {
    const url = '/database_individual';

    var data = {id: unique_id};

    const options = {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(data)
    };

    try {
        const response = await fetch(url, options);
        if (!response.ok) {
            throw new Error('Fehler bei der Ausführung der Abfrage: Status ' + response.status);
        }
        return await response.json();  // parsed the response body as JSON
    } catch (error) {
        throw new Error('Netzwerkfehler oder Fehler beim Parsen der Daten: ' + error.message);
    }
}


// parses, cleans  data and creates table
function createTable(tableData){
    // Generates the columns for Tabulator
   
    function generateColumns(vars) {
        // Add row number column
        let columns = [{
            title: "#",
            field: "_row",
            formatter: "rownum",
            width: 40,
            headerSort: false,
            download: false
        }];
    
        // Generate columns for the rest of the variables
        columns.push(...vars.map(varName => {
            let column = {
                title: varName.toUpperCase(),
                field: varName,
                sorter: "string",
                headerFilter: true,
            };
    
            // Add formatter to specific columns
            if (varName == "Name" || varName == "ID" || varName == "Unit" || varName == "Value") {
                column.formatter = function(cell) {
                    var value = cell.getValue();
                    if (value == "Download" || varName == "Unit") {
                        return "<span class='text-primary'><u>" + value + "</u></span>";
                    } else {
                        return value;
                    }
                };
            }
    
            return column;
        }));
    
        return columns;
    }

    function filterAndRemoveTuplesByUnit(tuplesList) {
    let filteredList = [];
    let remainingList = [];

    // Durchlaufen der Liste von Tupeln
    for (let i = 0; i < tuplesList.length; i++) {
        // Überprüfen, ob der erste Eintrag des Tupels auf "Unit" endet
        if (tuplesList[i][0].endsWith("Unit")) {
            filteredList.push(tuplesList[i]);  // Tupel in die gefilterte Liste packen
        } else {
            remainingList.push(tuplesList[i]);  // Tupel in die verbleibende Liste packen
        }
    }

    return { filteredList, remainingList };  // Beide Listen zurückgeben
    }

    function jsonToKeyValueTuples(jsonData) {
        let tupleList = [];
    
        // Durchlaufe alle Schlüssel im JSON-Objekt
        for (let key in jsonData) {
            if (jsonData.hasOwnProperty(key)) {
                // Erstelle ein Tupel [Schlüssel, Wert]
                let tuple = [key, jsonData[key]];
                tupleList.push(tuple);
            }
        }
        console.log(tupleList);
        return tupleList;
    }

    function appendUnitValues(filteredList, remainingList) {
        let updatedList = [...remainingList];  // Kopiere die verbleibende Liste
    
        // Hilfsfunktion, um Zahlen und Unterstriche zu entfernen
        function normalizeKey(key) {
            return key.replace(/[\d_]/g, "");  // Entfernt alle Unterstriche und Zahlen
        }
    
        // Durchlaufe alle Unit-Einträge
        for (let i = 0; i < filteredList.length; i++) {
            let unitKey = filteredList[i][0].replace("Unit", "");  // Entferne das "Unit"
            let normalizedUnitKey = normalizeKey(unitKey);  // Entferne Zahlen und Unterstriche
            let unitValue = filteredList[i][1];  // Der Wert der Unit
    
            // Durchlaufe die verbleibenden Einträge und prüfe, ob der Anfang des Schlüssels übereinstimmt
            for (let j = 0; j < updatedList.length; j++) {
                let normalizedRemainingKey = normalizeKey(updatedList[j][0]);  // Normalisiere den Schlüssel
    
                // Prüfen, ob die Schlüssel übereinstimmen, nachdem Zahlen und Unterstriche entfernt wurden
                if (normalizedRemainingKey.startsWith(normalizedUnitKey)) {
                    // Füge den Unit-Wert als drittes Element in das verbleibende Tupel ein
                    updatedList[j].push(unitValue);
                }
            }
        }
    
        return updatedList;
    }

    function formatTuples(triplesOrTuplesList, vars) {
        let formattedList = [];
    
        // Durchlaufe die Liste von Tupeln/Triples
        for (let i = 0; i < triplesOrTuplesList.length; i++) {
            let current = triplesOrTuplesList[i];
            
            // Überprüfen, ob es sich um ein Tupel oder Triple handelt
            if (current.length === 2) {
                // Falls es ein Tupel ist, füge eine leere Zeichenkette als drittes Element hinzu
                current.push("");
            }
    
            // Erstelle das Objekt mit den entsprechenden Variablen (var1, var2, var3)
            let formattedEntry = {
                [vars[0]]: current[0],
                [vars[1]]: current[1],
                [vars[2]]: current[2]
            };
    
            // Füge das formatierte Objekt zur Liste hinzu
            formattedList.push(formattedEntry);
        }
    
        return formattedList;
    }

    // Transforms JSON response data from Backend for Tabulator
    function transformData(tableData, vars) {
        let tupel_data = filterAndRemoveTuplesByUnit(jsonToKeyValueTuples(tableData));
        let jsondata = tupel_data.remainingList;
        let unitdata = tupel_data.filteredList;

        let data = appendUnitValues(unitdata, jsondata);
        console.log(data);
        return formatTuples(data, vars);
    }

    try {
        var columns = generateColumns(['Name', 'Value', 'Unit']);
        var data = transformData(tableData, ['Name', 'Value', 'Unit']);

        table.setColumns(columns);  // Setzt die dynamisch erzeugten Spalten
        table.setData(data);        // Setzt die transformierten Daten
    } catch (error) {
        document.getElementById('queryResults').textContent = 'Fehler beim Parsen der Daten: ' + error.message;
    }
}

// Functionality for downloading table
function downloadTable() {
    table.download("csv", "daten.csv");
}

var heading = document.getElementById('heading')
var searchField = document.getElementById('sparqlForm')
// Function for back button
function goBack() { 
    // Navigate back to the previous URL
    window.location.reload();
    // Hide the heading again after 1 sec
    setTimeout(function() {
        heading.classList.replace("d-md-flex", "d-none");
        searchField.classList.replace("d-none", "d-md-flex");
    }, 1000);
};

// Enable Browser back
window.onpopstate = function(event) {
    if (event.state) {
         create_query(event.state.view, event.state.name || "");
     }
};

// When the form is submitted
document.getElementById('sparqlForm').addEventListener('submit', function(e) {
    e.preventDefault();

    // Get the search term from the input field
    var search_term = document.getElementById('nameInput').value;

    // Define the URL and options for the fetch request
    const url = '/database_new';
    const options = {
        method: 'POST',
        headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
        body: 'nameInput=' + encodeURIComponent(search_term)
    };

    // Send a POST request to the server with the search term
    fetch(url, options)
        .then(response => {
            if (!response.ok) {
                throw new Error('Error executing the query: Status ' + response.status);
            }
            // Parse the response body as JSON and return it
            return response.json();
        })
        .then(data => {
            console.log(data)
            // Update the table with the new data
            table.replaceData(data); 
        })
        .catch(error => {
            console.error('Network error or error parsing the data: ' + error.message);
        });
});



