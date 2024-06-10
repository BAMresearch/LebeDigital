var table = null;

// Functionality for downloading table
function downloadTable() {
    table.download("csv", "daten.csv");
}

var heading = document.getElementById('heading')
// Function for back button
function goBack() { 
    history.back(); 
    // Hide the heading again after 1 sec
    setTimeout(function() {
        heading.classList.replace("d-md-flex", "d-none");
    }, 1000);
};

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

// Funktion zum Extrahieren der Werte einer spezifizierten Eigenschaft aus den Bindings
function extractValues(data, propertyName) {
    if (!data.message || !data.message.results || !data.message.results.bindings) {
        console.error("Fehler: Die Datenstruktur ist nicht wie erwartet.");
        return [];
    }

    // Extrahieren der 'value'-Werte für die spezifizierte Eigenschaft aus jedem Binding
    const propertyValues = data.message.results.bindings.map(binding => {
        if (binding[propertyName] && binding[propertyName].value) {
            return binding[propertyName].value; // Gibt den Wert zurück, wenn die Eigenschaft existiert
        } else {
            return null; // Gibt null zurück, wenn die Eigenschaft im Binding fehlt
        }
    });

    return propertyValues;
}

// Creates a Query based on the specific input from the user and executes it
async function create_query(selectedQueryType, enteredName){
    let query = null;

    // Show all Mixtures with Name and ID
    if (selectedQueryType === "Mischung" && enteredName === "") {
        query = `
        SELECT ?ID ?Name WHERE {
            ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/cpto/MaterialComposition>.
            ?s <http://purl.org/spar/datacite/hasIdentifier> ?o.
            ?s <http://purl.org/spar/datacite/hasIdentifier> ?p.
            ?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual>.
            ?o <https://w3id.org/pmd/co/value> ?Name.
            ?p <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual>.
            ?p <https://w3id.org/pmd/co/value> ?ID.
            FILTER(REGEX(str(?o), "/humanreadableID[^/]*$"))
            FILTER(REGEX(str(?p), "/ID[^/]*$"))
        }`;
        history.pushState({ view: "Mischung"}, "");
    }

    // Show searched Mixtures
    if (selectedQueryType === "Mischung" && enteredName !== ""){
        query = `SELECT ?ID ?Name WHERE {
            ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/cpto/MaterialComposition>.
            ?s <http://purl.org/spar/datacite/hasIdentifier> ?o.
            ?s <http://purl.org/spar/datacite/hasIdentifier> ?p.
            ?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual>.
            ?o <https://w3id.org/pmd/co/value> ?Name.
            ?p <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual>.
            ?p <https://w3id.org/pmd/co/value> ?ID.
            FILTER(REGEX(str(?o), "/humanreadableID[^/]*$"))
            FILTER(REGEX(str(?p), "/ID[^/]*$"))
            FILTER(CONTAINS(LCASE(?Name), "${enteredName}"))
            }`;
        history.pushState({ view: "Mischung",  name: enteredName }, "");
    }

    // Show all Info from a specific Mixture
    if (selectedQueryType === "Details" && enteredName !== ""){
        query = `
        SELECT ?Bestandteil ?Wert ?Einheit WHERE {
        ?g ?p "${enteredName}".
        BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
        ?Bestandteil <https://w3id.org/pmd/co/value> ?Wert.
        OPTIONAL { ?Bestandteil <https://w3id.org/pmd/co/unit> ?Einheit. }
        FILTER(STRENDS(STR(?Bestandteil), ?suffix))
        }`;
        history.pushState({ view: "Details", name: enteredName }, "");

        // Show the heading div after a delay so that data loading is complete
        setTimeout(function() {
            heading.classList.replace("d-none", "d-md-flex");
        }, 1000);

        var selectedName = document.getElementById('selectedName');
        selectedName.textContent = enteredName; 
    }

    // Show all Emodule with Name and ID
    if (selectedQueryType === "EModule" && enteredName === "") {
        query = `
        SELECT ?ID ?Mixture ?Name ?EModule WHERE {
            ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/pmd/co/ModulusOfElasticity>.
            ?o <https://w3id.org/pmd/co/input> ?s.
            ?o <http://purl.org/spar/datacite/hasIdentifier> ?p.
            ?o <http://purl.org/spar/datacite/hasIdentifier> ?d.
            ?o <http://purl.org/spar/datacite/hasIdentifier> ?m.
            FILTER(REGEX(str(?d), "/humanreadableID[^/]*$"))
            FILTER(REGEX(str(?p), "/ID[^/]*$"))
            FILTER(REGEX(str(?m), "/MixtureID[^/]*$"))
            OPTIONAL { ?p <https://w3id.org/pmd/co/value> ?ID. }
            OPTIONAL { ?m <https://w3id.org/pmd/co/value> ?Mixture. }
            OPTIONAL { ?d <https://w3id.org/pmd/co/value> ?Name. }
            OPTIONAL { ?s <https://w3id.org/pmd/co/value> ?EModule. }
        }`;
        history.pushState({ view: "EModule" }, "");
    }

    // Show all Info from a specific Emodule
    if (selectedQueryType === "EModule" && enteredName !== ""){
        query = `
        SELECT ?Bestandteil ?Wert ?Einheit WHERE {
        ?g ?p "${enteredName}".
        BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
        ?Bestandteil <https://w3id.org/pmd/co/value> ?Wert.
        OPTIONAL { ?Bestandteil <https://w3id.org/pmd/co/unit> ?Einheit. }
        FILTER(STRENDS(STR(?Bestandteil), ?suffix))
        }`;
        history.pushState({ view: "EModule", name: enteredName }, "");
    }

    // main logic for not extended queries
    try {
        const tableData = await executeSparqlQuery(query);
        createTable(tableData);
    } catch (error) {
        document.getElementById('queryResults').textContent = error.message;
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
            headerSort: false
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

// Define the function that you want to run
function runQuery(e) {
    if (e) e.preventDefault(); // Stops standard formular transmission

    // Get information from query type selector
    var queryTypeSelect = document.getElementById("queryType");
    var selectedQueryType = queryTypeSelect.value;
    console.log(selectedQueryType)

    // Get information from text input
    var nameInput = document.getElementById("nameInput");
    var enteredName = nameInput.value;
    console.log(enteredName)

    // creates a new table and adds placeholder
    table = new Tabulator("#resultsTable", {
    layout: "fitColumns",
    placeholder: "Loading Data...",
    pagination:"local",
    paginationSize:20,
    paginationSizeSelector:[10, 20, 50, 100],
    });


    // function for creating and executing Sparql query, based on input from the form
    create_query(selectedQueryType, enteredName)

    // Fügt den Event Listener für den Klick auf eine Zeile hinzu
    table.on("cellClick", function(e, cell) {
        // 'cell' ist das Zellen-Objekt, 'e' ist das Event-Objekt
        var value = cell.getValue(); // Holt den Wert der Zelle
        var field = cell.getField(); // Holt den Feldnamen der Spalte

        if (field == "Name" || field == "ID") { // Überprüft, ob die Spalte 'Name' ist
            document.getElementById('fileID').value = value;
            enteredName = value
            create_query("Details", enteredName);
        } else if (field == "Einheit") {
            var url = "https://qudt.org/vocab/unit/" + encodeURIComponent(value);
            window.open(url, "_blank");
        } else if (field == "Wert" && value == "Download") {
            // zuerst ID der mischung bekommen
            var id = document.getElementById('fileID').value
            // Sendet eine Anfrage an das Backend
            fetch(`/rawdownload?id=${encodeURIComponent(id)}`, {
                method: 'GET'
            }).then(response => {
                if (!response.ok) {
                    throw new Error('Network response was not ok');
                }
                return response.blob();
            }).then(blob => {
                // Erstellt einen Download-Link und klickt darauf
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
}

// Attach the function to the 'submit' event of the form
document.getElementById('sparqlForm').addEventListener('submit', runQuery);

// Also run the function when the page is loaded
window.onload = runQuery;

// Run the function when the dropdown value is changed
document.getElementById('queryType').addEventListener('change', runQuery);


// Enable Browser back
window.onpopstate = function(event) {
   if (event.state) {
        create_query(event.state.view, event.state.name || "");
    }
};

