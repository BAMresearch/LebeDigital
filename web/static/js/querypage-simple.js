var table = null;

// Functionality for downloading table
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
    // check if extended is checked
    const checkbox = document.getElementById('extendedCheckbox');
    let query = null;
    let response = null;
    let valuesList = null;

    if (checkbox.checked){
        // shows all information for every mixture
        if (selectedQueryType === "Mischung" && enteredName === "") {
            query = `
            SELECT ?ID WHERE {
              ?s <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <https://w3id.org/cpto/MaterialComposition>.
              ?s <http://purl.org/spar/datacite/hasIdentifier> ?p.
              ?p <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#NamedIndividual>.
              ?p <https://w3id.org/pmd/co/value> ?ID.
              FILTER(REGEX(str(?p), "/humanreadableID[^/]*$"))
            }`;
            response = await executeSparqlQuery(query);

            // extracts an array of every ID for every mixture
            valuesList = extractValues(response, "ID");
            console.log(valuesList);

            let allList = [];

            // goes through all mixtures in database
            for (const [index, value] of valuesList.entries()) {

                // looks up mixture
                query = `
                SELECT ?Bestandteil ?Wert WHERE {
                    ?g ?p "${value}".
                    BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
                    ?Bestandteil <https://w3id.org/pmd/co/value> ?Wert.
                    FILTER(STRENDS(STR(?Bestandteil), ?suffix))
                }`;

                response = await executeSparqlQuery(query);
                console.log(response)
                // extract all "Bestandteil"
                titleList = extractValues(response, "Bestandteil");

                // cleans data
                for (let i = 0; i < titleList.length; i++) {
                    titleList[i] = cleanEntry(titleList[i]);
                }

                // extract all "Bestandteil"
                valuesList = extractValues(response, "Wert");

                // for the first mixture extract column Names
                if (index === 0) {
                    // sets columnn name
                    console.log(titleList);
                    tableData = {message: {head: {vars: titleList}, results: {bindings: []}}};
                }

                let pushObject = {};
                for (let i = 0; i < titleList.length; i++) {
                    pushObject[titleList[i]] = { type: "literal", value: valuesList[i] };
                }

                console.log(pushObject);  // Dies zeigt das Objekt direkt in der Konsole
                allList.push(pushObject);  // Fügen Sie das Objekt zur Liste hinzu
            }
            tableData.message.results.bindings = allList;
            console.log(tableData);
            createTable(tableData);
        }

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
        if (selectedQueryType === "Mischung" && enteredName !== ""){
            query = `
            SELECT ?Bestandteil ?Wert ?Einheit WHERE {
            ?g ?p "${enteredName}".
            BIND(SUBSTR(STR(?g), STRLEN(STR(?g)) - 35) AS ?suffix)
            ?Bestandteil <https://w3id.org/pmd/co/value> ?Wert.
            OPTIONAL { ?Bestandteil <https://w3id.org/pmd/co/unit> ?Einheit. }
            FILTER(STRENDS(STR(?Bestandteil), ?suffix))
            }`;
        }

        // Show all Mixtures with Name and ID
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
        }

        // Show all Info from a specific Mixture
        if (selectedQueryType === "EModule" && enteredName !== ""){
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
    placeholder: "Daten werden geladen...",
    });


    // function for creating and executing Sparql query, based on input from the form
    create_query(selectedQueryType, enteredName)

    // Fügt den Event Listener für den Klick auf eine Zeile hinzu
    table.on("cellClick", function(e, cell) {
        // 'cell' ist das Zellen-Objekt, 'e' ist das Event-Objekt
        var value = cell.getValue(); // Holt den Wert der Zelle
        var field = cell.getField(); // Holt den Feldnamen der Spalte

        if (field == "Name" || field == "ID") { // Überprüft, ob die Spalte 'Name' ist
            document.getElementById('nameInput').value = value;  // Setzt den Wert ins Eingabefeld
            document.getElementById('submit').click(); // Automatisches Auslösen des Query-Buttons
        } else if (field == "Einheit") {
            var url = "https://qudt.org/vocab/unit/" + encodeURIComponent(value);
            window.open(url, "_blank");
        } else {
            console.log(`Geklickt auf Spalte: ${field} mit Wert: ${value}`);
        }
    });

});

