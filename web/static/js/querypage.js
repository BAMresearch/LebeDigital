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
            row[varName] = binding[varName] || "";
        });
        return row;
    });
}

function executeSparqlQuery(query, table) {
    var xhr = new XMLHttpRequest();
    const url = window.appConfig.urls.queryExecute;
    xhr.open('POST', url, true);
    xhr.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded');

    xhr.onload = function() {
        if (xhr.status >= 200 && xhr.status < 300) {
            try {
                var response = JSON.parse(xhr.responseText);
                console.log(response.message);
                var vars = Object.keys(response.message[0]);
                console.log(vars);

                var columns = generateColumns(vars);
                var data = transformData(response.message, vars);
                console.log(data);
                console.log(columns);
                table.setColumns(columns);  // Setzt die dynamisch erzeugten Spalten
                table.setData(data);        // Setzt die transformierten Daten

                document.getElementById('downloadBtn').classList.replace("d-none", "d-md-flex"); //show download button
            } catch (error) {
                console.log(error.message)
                //document.getElementById('queryResults').innerHTML = 'Please insert a valid query. ' + error.message;
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
    runQuery()
});

// Function to run the query
function runQuery() {
    var query = document.getElementById('sparqlQuery').value;

    table = new Tabulator("#resultsTable", {
        layout: "fitColumns",
        placeholder: "Loading Data...",
        pagination: "local",
        paginationSize: 10,
        paginationSizeSelector: [10, 50, 100],
    });

    // Execute the query using the defined function
    executeSparqlQuery(query, table);
}

// Event-Handler für den Button, um den vorgefertigten Text einzufügen
document.getElementById('insertTextBtn').addEventListener('click', function() {
    var predefinedText = 'SELECT ?s ?p ?o WHERE { ?s ?p ?o }';
    document.getElementById('sparqlQuery').value = predefinedText;
});

// Define the SPARQL queries
const queries = {
    query1: `SELECT ?humanReadableID ?waterCementRatio
        WHERE {
        ?mixture a <https://w3id.org/cpto/MaterialComposition> .
        ?mixture <http://purl.org/spar/datacite/hasIdentifier> ?idNode .
        ?idNode a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                <https://w3id.org/pmd/co/value> ?humanReadableID .
        FILTER(CONTAINS(STR(?idNode), "humanreadableID"))

        ?mixture <https://w3id.org/pmd/co/characteristic> ?wcrNode .
        ?wcrNode a <https://w3id.org/cpto/WaterCementRatio> ;
                <https://w3id.org/pmd/co/value> ?waterCementRatio .
        }`,
    query2: `SELECT ?humanReadableID ?cementContent ?cementType
        WHERE {
        ?mixture a <https://w3id.org/cpto/MaterialComposition> .
        ?mixture <http://purl.org/spar/datacite/hasIdentifier> ?idNode .
        ?idNode a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                <https://w3id.org/pmd/co/value> ?humanReadableID .
        FILTER(CONTAINS(STR(?idNode), "humanreadableID"))

        ?mixture <https://w3id.org/pmd/co/characteristic> ?materialComp .
        ?materialComp <https://w3id.org/pmd/co/composedOf> ?cement .
        ?cement <https://w3id.org/pmd/co/composedOf> ?cementTypeNode .
        ?cementTypeNode a <https://w3id.org/cpto/Cement> ;
                        <https://w3id.org/pmd/co/value> ?cementType .

        ?cement <https://w3id.org/pmd/co/characteristic> ?contentNode .
        ?contentNode a <https://w3id.org/cpto/Content> ;
                    <https://w3id.org/pmd/co/value> ?cementContent .
        }`,
    query3: `SELECT ?humanReadableID ?admixtureName ?admixtureContent
        WHERE {
        ?mixture a <https://w3id.org/cpto/MaterialComposition> .
        ?mixture <http://purl.org/spar/datacite/hasIdentifier> ?idNode .
        ?idNode a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                <https://w3id.org/pmd/co/value> ?humanReadableID .
        FILTER(CONTAINS(STR(?idNode), "humanreadableID"))

        ?mixture <https://w3id.org/pmd/co/characteristic> ?materialComp .
        ?materialComp <https://w3id.org/pmd/co/composedOf> ?admixture .
        ?admixture <https://w3id.org/pmd/co/composedOf> ?admixtureTypeNode .
        ?admixtureTypeNode a <https://w3id.org/cpto/Admixture> ;
                            <https://w3id.org/pmd/co/value> ?admixtureName .

        ?admixture <https://w3id.org/pmd/co/characteristic> ?contentNode .
        ?contentNode a <https://w3id.org/cpto/Content> ;
                    <https://w3id.org/pmd/co/value> ?admixtureContent .
        }`,
    query4: `SELECT ?humanReadableID ?compressiveStrength ?WaterCementRatio
        WHERE {
        ?specimen a <https://w3id.org/pmd/co/Specimen> .
        ?specimen <http://purl.org/spar/datacite/hasIdentifier> ?idNode .
        ?idNode a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                <https://w3id.org/pmd/co/value> ?humanReadableID .
        FILTER(CONTAINS(STR(?idNode), "humanreadableID"))

        OPTIONAL {
            ?specimen <https://w3id.org/pmd/co/input> ?csNode .
            ?csNode a <https://w3id.org/cpto/ConcreteCompressiveStrength> ;
                    <https://w3id.org/pmd/co/value> ?compressiveStrength .
        }

        OPTIONAL {
            ?MaterialComposition a <https://w3id.org/cpto/MaterialComposition> .
			?characteristic a <https://w3id.org/cpto/WaterCementRatio> .
			?MaterialComposition <https://w3id.org/pmd/co/characteristic> ?characteristic .
			?characteristic <https://w3id.org/pmd/co/value> ?WaterCementRatio  .
		}
        }`,
    query5: `SELECT ?humanReadableID ?compressiveStrength ?elasticModulus
        WHERE {
        ?specimen a <https://w3id.org/pmd/co/Specimen> .
        ?specimen <http://purl.org/spar/datacite/hasIdentifier> ?idNode .
        ?idNode a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                <https://w3id.org/pmd/co/value> ?humanReadableID .
        FILTER(CONTAINS(STR(?idNode), "humanreadableID"))

        OPTIONAL {
            ?specimen <https://w3id.org/pmd/co/input> ?csNode .
            ?csNode a <https://w3id.org/cpto/ConcreteCompressiveStrength> ;
                    <https://w3id.org/pmd/co/value> ?compressiveStrength .
        }

        OPTIONAL {
            ?specimen <https://w3id.org/pmd/co/input> ?emNode .
            ?emNode a <https://w3id.org/pmd/co/ModulusOfElasticity> ;
                    <https://w3id.org/pmd/co/value> ?elasticModulus .
        }
        }`,
    query6: `SELECT ?humanReadableID ?E_Module ?WaterCementRatio
        WHERE {
        ?specimen a <https://w3id.org/pmd/co/Specimen> .
        ?specimen <http://purl.org/spar/datacite/hasIdentifier> ?idNode .
        ?idNode a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                <https://w3id.org/pmd/co/value> ?humanReadableID .
        FILTER(CONTAINS(STR(?idNode), "humanreadableID"))

        OPTIONAL {
            ?specimen <https://w3id.org/pmd/co/input> ?csNode .
            ?csNode a <https://w3id.org/pmd/co/ModulusOfElasticity> ;
                    <https://w3id.org/pmd/co/value> ?E_Module .
        }

        OPTIONAL {
            ?MaterialComposition a <https://w3id.org/cpto/MaterialComposition> .
			?characteristic a <https://w3id.org/cpto/WaterCementRatio> .
			?MaterialComposition <https://w3id.org/pmd/co/characteristic> ?characteristic .
			?characteristic <https://w3id.org/pmd/co/value> ?WaterCementRatio  .
		}
        }`
};
// Event listener for all query items
document.querySelectorAll('.query-item').forEach(function(item) {
    item.addEventListener('click', function(e) {
        e.preventDefault();  // Prevent default action of anchor tag

        // Remove 'active' class from all query items
        document.querySelectorAll('.query-item').forEach(function(el) {
            el.classList.remove('active');
        });

        // Add 'active' class to the clicked item
        e.currentTarget.classList.add('active');

        // Get the query key from data attribute and set the query in textarea
        const queryKey = e.currentTarget.getAttribute('data-query');
        document.getElementById('sparqlQuery').value = queries[queryKey];

        // Automatically run the query after inserting it
        runQuery();
    });
});



