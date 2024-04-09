var table = null;

// Functionality for downloading table
function downloadTable() {
    table.download("csv", "daten.csv");
}

// Creates a Query based on the specific input from the user and executes it
async function create_query(selectedQueryType, enteredName){
    // check if extended is checked
    const checkbox = document.getElementById('extendedCheckbox');
    let query = null;

    if (checkbox.checked){
        // code
    } else {
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
        if (selectedQueryType === "Mischung" && enteredName !== "" && !checkbox.checked){
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

}


// parses, cleans  data and creates table
function createTable(tableData){

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
                // Prüft, ob das Objekt existiert, bevor auf .value zugegriffen wird
                let rawValue = binding[varName] ? binding[varName].value : "";
                // Reinigt den Wert vor der Zuweisung zum row-Objekt
                let cleanedValue = cleanEntry(rawValue);
                row[varName] = cleanedValue;
            });
            return row;
        });
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


// Form event Listener
document.getElementById('sparqlForm').addEventListener('submit', function(e) {
    // Starts when "Run Query" is pressed

    e.preventDefault(); // Stops standard formular transmission

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
    placeholder: "Daten werden geladen..."
    });

    // function for creating and executing Sparql query, based on input from the form
    create_query(selectedQueryType, enteredName)

});