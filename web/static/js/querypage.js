var table = null;
var myChart = null;

function initializeOrUpdateTable(data, vars) {
    const columns = generateColumns(vars);

    if (!table) {
        // Initialize the table with client-side pagination
        table = new Tabulator("#resultsTable", {
            columns: columns,
            layout: "fitColumns",
            placeholder: "Loading Data...",
            pagination: true,           // Enable client-side pagination
            paginationSize: 10,         // Number of rows per page
            paginationInitialPage: 1,   // Start on the first page
            data: data.message,         // Load all data at once
        });
    } else {
        // Update columns and data if table already exists
        table.setColumns(columns);
        table.setData(data.message);
    }
}

// Helper function to fetch data asynchronously
async function executeSparqlQuery(query) {
    const url = window.appConfig.urls.queryExecute;
    try {
        const response = await fetch(url, {
            method: 'POST',
            headers: { 'Content-Type': 'application/x-www-form-urlencoded' },
            body: `query=${encodeURIComponent(query)}`
        });

        if (!response.ok) throw new Error('Network response was not ok');

        const result = await response.json();
        console.log(result)
        if (!result.message || result.message.length === 0) {
            document.getElementById('queryResults').innerHTML = 'No results found';
            clearTableAndChart();
            return null;
        }

        return result;

    } catch (error) {
        console.error('Error fetching data:', error);
        document.getElementById('queryResults').innerHTML = 'Error fetching data';
        clearTableAndChart();
        return null;
    }
}

// Clear chart and table on error or no data
function clearTableAndChart() {
    clearTable();
    clearChart();
}

function clearChart() {
    if (myChart) {
        myChart.destroy();
        myChart = null;
    }
}

function clearTable(){
    if (table) table.clearData();
}


// Generate dynamic columns based on query variables
function generateColumns(vars) {
    let columns = [{
        title: "#",
        formatter: function(cell) {
            let table = cell.getTable();
            let row = cell.getRow();
            let page = table.getPage();
            let pageSize = table.getPageSize();
            return ((page - 1) * pageSize) + row.getPosition(true) + 1;
        },
        width: 40,
        headerSort: false,
        download: false
    }];
    
    const dataColumns = vars.map(varName => ({
        title: varName.toUpperCase(),
        field: varName,
        sorter: "string",
        headerFilter: true
    }));

    return columns.concat(dataColumns);
}


function updateChart(data, vars) {
    const ctx = document.getElementById('myChart').getContext('2d');
    const queryType = determineQueryType(vars);

    if (!queryType || !queryVisualizationSettings[queryType]) {
        clearChart(); // Ensure the chart is cleared if the query type is not valid
        return;
    }

    const chartConfig = queryVisualizationSettings[queryType].config(data);

    if (myChart) {
        // Destroy the previous chart instance if the type has changed
        if (myChart.config.type !== chartConfig.type) {
            clearChart();
        } else {
            // Update chart data and options if chart already exists
            myChart.data = chartConfig.data;
            myChart.options = chartConfig.options;
            myChart.update();
            return;
        }
    }

    // Create a new chart instance
    myChart = new Chart(ctx, chartConfig);
}

// Main function to run the query and update UI
async function runQuery() {
    const query = document.getElementById('sparqlQuery').value;
    const data = await executeSparqlQuery(query, 1, 10); // Fetch first page with page size 10

    if (!data) return;

    const vars = Object.keys(data.message[0]);
    initializeOrUpdateTable(data, vars, data.totalResults);  // Initialize or update table
    updateChart(data.message, vars);

    document.getElementById('downloadBtn').classList.replace("d-none", "d-md-flex");
    document.getElementById('queryResults').innerHTML = '';
}

// Event listener for the form submission
document.getElementById('sparqlForm').addEventListener('submit', function (e) {
    e.preventDefault();
    runQuery();
});

// Button to insert predefined text into the query input
document.getElementById('insertTextBtn').addEventListener('click', function () {
    const predefinedText = 'SELECT ?s ?p ?o WHERE { ?s ?p ?o }';
    document.getElementById('sparqlQuery').value = predefinedText;
});

// Download table as CSV
function downloadTable() {
    if (table) {
        table.download("csv", "daten.csv");
    }
}

// Determine query type based on variables present in data
function determineQueryType(vars) {
    console.log(vars)
    if (vars.includes('waterCementRatio') && vars.length === 2) {
        return 'waterCementRatio';
    } else if (vars.includes('WaterCementRatio') && vars.includes('CompressiveStrength')) {
        return 'comst-wc';
    } else if (vars.includes('WaterCementRatio') && vars.includes('E_Module')) {
        return 'emod-wc';
    } else if (vars.includes('InputCompressiveStrength') && vars.includes('elasticModulus')) {
        return 'emod-comst';
    } else if (vars.includes('admixtureName')) {
        return 'admixture';
    } else if (vars.includes('cementType')) {
        return 'cement';
    }
    return null;
}

// Visualization settings for different query types
const queryVisualizationSettings = {
    // Define visualization configurations here, similar to before
    waterCementRatio: {
        type: 'bar',
        config: (data) => ({
            type: 'bar',
            data: {
                labels: data.map(d => d.humanReadableID),
                datasets: [{
                    label: 'Water-Cement Ratio',
                    data: data.map(d => parseFloat(d.waterCementRatio)),
                    backgroundColor: data.map(d => `hsl(${Math.random() * 360}, 70%, 50%)`),
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true
            }
        })
    },
    'admixture': {
        type: 'bar',
        config: (data) => ({
            type: 'bar',
            data: {
                labels: data.map(d => d.humanReadableID),
                datasets: [{
                    label: 'Admixture Content',
                    data: data.map(d => parseFloat(d.admixtureDensity)),
                    backgroundColor: data.map(d => `hsl(${Math.random() * 360}, 70%, 50%)`),
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const dataPoint = data[context.dataIndex];
                                return `${dataPoint.admixtureName}: ${dataPoint.admixtureDensity}`;
                            }
                        }
                    }
                },
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Density'
                        }
                    }
                }
            }
        })
    },

    'cement': {
        type: 'bar',
        config: (data) => ({
            type: 'bar',
            data: {
                labels: data.map(d => d.humanReadableID),
                datasets: [{
                    label: 'Cement Content',
                    data: data.map(d => parseFloat(d.cementContent)),
                    backgroundColor: data.map(d => `hsl(${Math.random() * 360}, 70%, 50%)`),
                    borderWidth: 1
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const dataPoint = data[context.dataIndex];
                                return `${dataPoint.cementType}: ${dataPoint.cementContent}`;
                            }
                        }
                    }
                },
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Content'
                        }
                    }
                }
            }
        })
    },
    'comst-wc': {
        type: 'scatter',
        config: (data) => ({
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Water-Cement Ratio vs Compressive Strength',
                    data: data.map(d => ({
                        x: parseFloat(d.WaterCementRatio),
                        y: parseFloat(d.CompressiveStrength),
                        humanReadableID: d.humanReadableID  // Store ID for tooltip
                    })),
                    borderColor: 'rgb(75, 192, 192)',
                    backgroundColor: 'rgba(75, 192, 192, 0.7)',
                    pointRadius: 8,
                    showLine: false
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const dataPoint = context.raw;
                                return [
                                    `ID: ${dataPoint.humanReadableID}`,
                                    `WC Ratio: ${dataPoint.x}`,
                                    `Strength: ${dataPoint.y} MPa`
                                ];
                            }
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'Water-Cement Ratio'
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'Compressive Strength (MPa)'
                        }
                    }
                }
            }
        })
    },
    'emod-wc': {
        type: 'scatter',
        config: (data) => ({
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Water-Cement Ratio vs E-Module',
                    data: data.map(d => ({
                        x: parseFloat(d.WaterCementRatio),
                        y: parseFloat(d.E_Module),
                        humanReadableID: d.humanReadableID  // Store ID for tooltip
                    })),
                    borderColor: 'rgb(54, 162, 235)',
                    backgroundColor: 'rgba(54, 162, 235, 0.7)',
                    pointRadius: 6,
                    showLine: false
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const dataPoint = context.raw;
                                return [
                                    `ID: ${dataPoint.humanReadableID}`,
                                    `WC Ratio: ${dataPoint.x}`,
                                    `Strength: ${dataPoint.y} N/mm²`
                                ];
                            }
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'Water-Cement Ratio'
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'E-Module (N/mm²)'
                        }
                    }
                }
            }
        })
    },
    'emod-comst': {
        type: 'scatter',
        config: (data) => ({
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Input Compressive Strength vs E-Module',
                    data: data.map(d => ({
                        x: parseFloat(d.InputCompressiveStrength),
                        y: parseFloat(d.elasticModulus),
                        humanReadableID: d.humanReadableID  // Store ID for tooltip
                    })),
                    borderColor: 'rgb(255, 99, 132)',
                    backgroundColor: 'rgba(255, 99, 132, 0.7)',
                    pointRadius: 6,
                    showLine: false
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const dataPoint = context.raw;
                                return [
                                    `ID: ${dataPoint.humanReadableID}`,
                                    `WC Ratio: ${dataPoint.x}`,
                                    `Strength: ${dataPoint.y} N/mm²`
                                ];
                            }
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'Input Compressive Strength'
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'E-Module (N/mm²)'
                        }
                    }
                }
            }
        })
    }
};

// Event listener for running predefined queries
document.querySelectorAll('.query-item').forEach(function (item) {
    item.addEventListener('click', function (e) {
        e.preventDefault();
        document.querySelectorAll('.query-item').forEach(el => el.classList.remove('active'));
        e.currentTarget.classList.add('active');
        const queryKey = e.currentTarget.getAttribute('data-query');
        document.getElementById('sparqlQuery').value = queries[queryKey];
        runQuery();
    });
});

// Define the SPARQL queries (sample queries)
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
    query2:`SELECT ?humanReadableID ?cementContent ?cementType
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
    query3:  `SELECT ?humanReadableID ?admixtureName ?admixtureDensity
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

        ?admixture <https://w3id.org/pmd/co/characteristic> ?densityNode .
        ?densityNode a <https://w3id.org/cpto/RelativeDensity> ;
                <https://w3id.org/pmd/co/value> ?admixtureDensity .
        }`,
    query4: `SELECT DISTINCT ?humanReadableID ?WaterCementRatio ?CompressiveStrength WHERE { 
            # Material Composition and WaterCementRatio
            ?MaterialComposition a <https://w3id.org/cpto/MaterialComposition> ;
                                <http://purl.org/spar/datacite/hasIdentifier>/<https://w3id.org/pmd/co/value> ?ID_1 ;
                                <https://w3id.org/pmd/co/characteristic> ?wcrNode .
            
            ?wcrNode <https://w3id.org/pmd/co/value> ?WaterCementRatio .
            
            # Ensure WCR type
            FILTER EXISTS {
                ?wcrNode a <https://w3id.org/cpto/WaterCementRatio>
            }
            
            # Specimen and measurements
            ?Specimen a <https://w3id.org/pmd/co/Specimen> ;
                    <http://purl.org/spar/datacite/hasIdentifier> ?hasIdentifier_2 .
            
            # Connect specimen to material composition
            ?Specimen <http://purl.org/spar/datacite/hasIdentifier> ?providedId .
            ?providedId a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                        <https://w3id.org/pmd/co/value> ?ID_1 .
            
            # Get humanReadableID and Compressive Strength
            ?hasIdentifier_2 <https://w3id.org/pmd/co/value> ?humanReadableID .
            
            ?ComSt a <https://w3id.org/cpto/ConcreteCompressiveStrength> ;
                <http://purl.org/spar/datacite/hasIdentifier> ?hasIdentifier_2 ;
                <https://w3id.org/pmd/co/value> ?CompressiveStrength .
        }`,
    query5: `SELECT DISTINCT ?humanReadableID ?InputCompressiveStrength ?elasticModulus WHERE {
        # Retrieve the specimen and its associated human-readable ID
        ?specimen a <https://w3id.org/pmd/co/Specimen> .
        ?specimen <http://purl.org/spar/datacite/hasIdentifier> ?idNode .
        ?idNode a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                <https://w3id.org/pmd/co/value> ?humanReadableID .
        FILTER(CONTAINS(STR(?idNode), "humanreadableID"))

        # Retrieve Modulus of Elasticity (elasticModulus) linked as an input to the specimen
        ?specimen <https://w3id.org/pmd/co/input> ?elasticityNode .
        ?elasticityNode a <https://w3id.org/pmd/co/ModulusOfElasticity> ;
                        <https://w3id.org/pmd/co/value> ?elasticModulus .

        # Link specimen to its associated experiment info with Input Compressive Strength
        ?experimentInfo a <https://w3id.org/cpto/DeterminationOfSecantModulusOfElasticity> ;
                        <https://w3id.org/pmd/co/output> ?outputNode .
        ?outputNode a <https://w3id.org/pmd/co/InputCompressiveStrength> ;
                    <https://w3id.org/pmd/co/value> ?InputCompressiveStrength .

        # Ensure the experiment info is directly associated with the specimen
        ?elasticityNode <https://w3id.org/pmd/co/output> ?experimentInfo .
        }`,
    query6: `SELECT DISTINCT ?humanReadableID ?WaterCementRatio ?E_Module WHERE { 
            # Material Composition and WaterCementRatio
            ?MaterialComposition a <https://w3id.org/cpto/MaterialComposition> ;
                                <http://purl.org/spar/datacite/hasIdentifier>/<https://w3id.org/pmd/co/value> ?ID_1 ;
                                <https://w3id.org/pmd/co/characteristic> ?wcrNode .
            
            ?wcrNode <https://w3id.org/pmd/co/value> ?WaterCementRatio .
            
            # Ensure WCR type
            FILTER EXISTS {
                ?wcrNode a <https://w3id.org/cpto/WaterCementRatio>
            }
            
            # Specimen and measurements
            ?Specimen a <https://w3id.org/pmd/co/Specimen> ;
                    <http://purl.org/spar/datacite/hasIdentifier> ?hasIdentifier_2 .
            
            # Connect specimen to material composition
            ?Specimen <http://purl.org/spar/datacite/hasIdentifier> ?providedId .
            ?providedId a <https://w3id.org/pmd/co/ProvidedIdentifier> ;
                        <https://w3id.org/pmd/co/value> ?ID_1 .
            
            # Get humanReadableID and E-Module
            ?hasIdentifier_2 <https://w3id.org/pmd/co/value> ?humanReadableID .
            
            ?ElasticM a <https://w3id.org/pmd/co/ModulusOfElasticity> ;
                    <http://purl.org/spar/datacite/hasIdentifier> ?hasIdentifier_2 ;
                    <https://w3id.org/pmd/co/value> ?E_Module .
        }`,
};
