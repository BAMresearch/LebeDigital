var table = null;
var query = null;

// Hinzufügen eines Buttons zum Herunterladen als CSV
function downloadTable() {
    table.download("csv", "daten.csv");
}

function create_query(selectedQueryType, enteredName){
    if (selectedQueryType === "Mischung" && enteredName === "") {
        query = `
        SELECT ?b WHERE {
          ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/cpto/MaterialComposition>.
          ?s <http://purl.org/spar/datacite/hasIdentifier> ?o.
          ?o <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual>.
          ?o <https://w3id.org/pmd/co/value> ?b.
          FILTER(REGEX(str(?o), "/humanreadableID[^/]*$"))
        }`;
        console.log("a")
    }
    console.log("a")
    if (selectedQueryType === "Mischung" && enteredName !== ""){
        query = `
        SELECT ?a ?b WHERE {
        ?g ?p "${enteredName}".
        BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
        ?a <https://w3id.org/pmd/co/value> ?b.
        FILTER(STRENDS(STR(?a), ?suffix))
        }`;
        console.log(query)
    }
}
function generateColumns(vars) {
    return vars.map(varName => ({
        title: varName.toUpperCase(),  // Spaltentitel als Großbuchstaben der Variablennamen
        field: varName,
        sorter: "string",
        headerFilter: true  // Optional: Fügt Filtermöglichkeiten zu jeder Spalte hinzu
    }));
}

function transformData(bindings, vars) {
    return bindings.map(binding => {
        let row = {};
        vars.forEach(varName => {
            // Extrahiert den Wert für jede Variable in jeder Bindung
            row[varName] = binding[varName].value;
        });
        return row;
    });
}

function executeSparqlQuery() {
    var xhr = new XMLHttpRequest();
    xhr.open('POST', '/queryexec', true);
    xhr.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded');

    xhr.onload = function() {
        if (xhr.status >= 200 && xhr.status < 300) {
            try {
                var response = JSON.parse(xhr.responseText);
                var vars = response.message.head.vars;
                var bindings = response.message.results.bindings;

                var columns = generateColumns(vars);
                var data = transformData(bindings, vars);

                table.setColumns(columns);  // Setzt die dynamisch erzeugten Spalten
                table.setData(data);        // Setzt die transformierten Daten
            } catch (error) {
                document.getElementById('queryResults').innerHTML = 'Fehler beim Parsen der Daten: ' + error.message;
            }
        } else {
            document.getElementById('queryResults').innerHTML = 'Fehler bei der Ausführung der Abfrage: Status ' + xhr.status;
        }
    };

    xhr.onerror = function() {
        document.getElementById('queryResults').innerHTML = 'Netzwerkfehler bei der Anfrage.';
    };

    xhr.send('query=' + encodeURIComponent(query));
}



document.getElementById('sparqlForm').addEventListener('submit', function(e) {
    e.preventDefault(); // Standard-Formular-Submit unterbrechen
    var queryTypeSelect = document.getElementById("queryType");
    var selectedQueryType = queryTypeSelect.value;
    console.log(selectedQueryType)

    var nameInput = document.getElementById("nameInput");
    var enteredName = nameInput.value;

    console.log(enteredName)

    create_query(selectedQueryType, enteredName)

    table = new Tabulator("#resultsTable", {
    layout: "fitColumns",
    placeholder: "Daten werden geladen..."
    });

    // Führen Sie die Abfrage mit der definierten Funktion aus
    executeSparqlQuery();
});