var table = null;

// Hinzufügen eines Buttons zum Herunterladen als CSV
function downloadTable() {
    table.download("csv", "daten.csv");
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
            row[varName] = binding[varName] ? binding[varName].value : "";
        });
        return row;
    });
}

function executeSparqlQuery(query, table) {
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

                document.getElementById('downloadBtn').style.display = 'block'; //show download button
            } catch (error) {
                document.getElementById('queryResults').innerHTML = 'Please insert a valid query. ' + error.message;
            }
        } else {
            document.getElementById('queryResults').innerHTML = 'Failed to load data: Status ' + xhr.status;
        }
    };

    xhr.onerror = function() {
        document.getElementById('queryResults').innerHTML = 'Network Error';
    };

    xhr.send('query=' + encodeURIComponent(query));
}



document.getElementById('sparqlForm').addEventListener('submit', function(e) {
    e.preventDefault(); // Standard-Formular-Submit unterbrechen
    var query = document.getElementById('sparqlQuery').value;

    table = new Tabulator("#resultsTable", {
    layout: "fitColumns",
    placeholder: "Daten werden geladen..."
    });

    // Führen Sie die Abfrage mit der definierten Funktion aus
    executeSparqlQuery(query, table);
});

// Event-Handler für den Button, um den vorgefertigten Text einzufügen
document.getElementById('insertTextBtn').addEventListener('click', function() {
    var predefinedText = 'SELECT ?s ?p ?o WHERE { ?s ?p ?o }';
    document.getElementById('sparqlQuery').value = predefinedText;
});