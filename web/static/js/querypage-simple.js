var table = null;
var query = null;

// Hinzufügen eines Buttons zum Herunterladen als CSV
function downloadTable() {
    table.download("csv", "daten.csv");
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


// Creates a Query based on the specific input from the user
function create_query(selectedQueryType, enteredName){
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
    }

    // Show all Info from a specific Mixture
    if (selectedQueryType === "Mischung" && enteredName !== ""){
        query = `
        SELECT ?Bestandteil ?Wert WHERE {
        ?g ?p "${enteredName}".
        BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
        ?Bestandteil <https://w3id.org/pmd/co/value> ?Wert.
        FILTER(STRENDS(STR(?Bestandteil), ?suffix))
        }`;
    }
}

// Generates the columns for Tabulator
function generateColumns(vars) {
    return vars.map(varName => ({
        title: varName.toUpperCase(),  // Spaltentitel als Großbuchstaben der Variablennamen
        field: varName,
        sorter: "string",
        headerFilter: true  // Optional: Fügt Filtermöglichkeiten zu jeder Spalte hinzu
    }));
}

// Transforms JSON response data from Backend for Tabulator
function transformData(bindings, vars) {
    return bindings.map(binding => {
        let row = {};
        vars.forEach(varName => {
            // Extrahiert den Wert für jede Variable in jeder Bindung
            let rawValue = binding[varName].value;
            // Reinigt den Wert vor der Zuweisung zum row-Objekt
            let cleanedValue = cleanEntry(rawValue);
            row[varName] = cleanedValue;
        });
        return row;
    });
}

// Send Sparql Query to backend and transform results
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


// Form event Listener
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