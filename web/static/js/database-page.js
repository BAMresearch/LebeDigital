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

    if (field == "filename" || field == "Unique_ID" || field == "type") {
        rowData = cell.getRow().getData();
        document.getElementById('fileID').value = rowData["Unique_ID"];
        var enteredName = rowData["filename"];
        var enteredType = rowData["type"];

        // Remove file extension from filename
        //var filenameWithoutExtension = enteredName.split('.').slice(0, -1).join('.');  
        
        if(enteredType !== "Mixture"){
            create_query(rowData["Unique_ID"], enteredType);
        }
        else{
            create_query(enteredName, enteredType);
        }

        // Show the heading div after a delay so that data loading is complete
        setTimeout(function() {
            heading.classList.replace("d-none", "d-md-flex");
            searchField.classList.replace("d-md-flex", "d-none");
        }, 1000);

        var selectedName = document.getElementById('selectedName');
        selectedName.textContent = filenameWithoutExtension; 
        
    } else if (field == "Einheit") {
        var url = "https://qudt.org/vocab/unit/" + encodeURIComponent(value);
        window.open(url, "_blank");
    } else if (field == "Wert" && value == "Download") {
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
        console.log(`Geklickt auf Spalte: ${field} mit Wert: ${value}`);
    }
});

// Creates a Query based on the specific input from the user and executes it
async function create_query(enteredName, fileType){
    let query = null;

    // Show all Info from a specific Mixture
    if (fileType === "Mixture" && enteredName !== ""){
        query = `
        SELECT ?Bestandteil ?Wert ?Einheit WHERE {
        ?g ?p "${enteredName}".
        BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
        ?Bestandteil <https://w3id.org/pmd/co/value> ?Wert.
        OPTIONAL { ?Bestandteil <https://w3id.org/pmd/co/unit> ?Einheit. }
        FILTER(STRENDS(STR(?Bestandteil), ?suffix))
        }`;
        
    }

    // Show all Info from a specific Emodule
    if (fileType === "EModule" && enteredName !== ""){
        query = `
        SELECT ?Bestandteil ?Wert ?Einheit WHERE {
        ?g ?p "${enteredName}".
        BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
        ?Bestandteil <https://w3id.org/pmd/co/value> ?Wert.
        OPTIONAL { ?Bestandteil <https://w3id.org/pmd/co/unit> ?Einheit. }
        FILTER(STRENDS(STR(?Bestandteil), ?suffix))
        }`;
    }

    // main logic for not extended queries
    try {
        const tableData = await executeSparqlQuery(query);
        createTable(tableData);
    } catch (error) {
        document.getElementById('queryResults').textContent = error.message;
    }
}

// executes the sparql query and returns result as a promise
async function executeSparqlQuery(query) {
    const url = '/queryexec';
    const options = {
        method: 'POST',
        headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
        body: 'query=' + encodeURIComponent(query)
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
            if (varName == "Name" || varName == "ID" || varName == "Einheit" || varName == "Wert") {
                column.formatter = function(cell) {
                    var value = cell.getValue();
                    if (value == "Download" || varName == "Einheit") {
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
    
    // Transforms JSON response data from Backend for Tabulator
    function transformData(bindings, vars) {
        return bindings.map(binding => {
            let row = {};
            vars.forEach(varName => {
                // Prüft, ob das Objekt existiert, bevor auf .value zugegriffen wird
                let rawValue = binding[varName] ? binding[varName].value : "";
                // Reinigt den Wert vor der Zuweisung zum row-Objekt
                let cleanedValue = cleanEntry(rawValue);
                row[varName] = cleanedValue;
            });
            return row;
        });
    }

    try {
        var vars = tableData.message.head.vars;
        var bindings = tableData.message.results.bindings;

        var columns = generateColumns(vars);
        var data = transformData(bindings, vars);
        table.setColumns(columns);  // Setzt die dynamisch erzeugten Spalten
        table.setData(data);        // Setzt die transformierten Daten
    } catch (error) {
        document.getElementById('queryResults').textContent = 'Fehler beim Parsen der Daten: ' + error.message;
    }
}

// Cleans the data retrieved from Sparql Query (Removes URI and UUID4)
function cleanEntry(entry) {
    // Überprüfung auf UUID4 am Ende und Entfernen
    const uuidRegex = /_[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$/;
    if (uuidRegex.test(entry)) {
        entry = entry.replace(uuidRegex, '');
    }

    // Überprüfung, ob der Eintrag mit einer URL beginnt und Entfernen bis zum letzten "/"
    const urlRegex = /^(https?:\/\/[^\/]+\/.*)$/;
    if (urlRegex.test(entry)) {
        const lastSlashIndex = entry.lastIndexOf('/');
        if (lastSlashIndex > -1) {
            entry = entry.substring(lastSlashIndex + 1);
        }
    }

    return entry;
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



