var table = null;

// Hinzufügen eines Buttons zum Herunterladen als CSV
function downloadTable() {
    table.download("csv", "daten.csv");
}

// Helper function to extract variable order from SPARQL query
function extractVariableOrder(query) {
    try {
        // Look for SELECT clause and capture variables
        const selectMatch = query.match(/SELECT\s+([?]\w+\s+)+/i);
        if (!selectMatch) return null;
        
        // Extract variables (words starting with ?)
        const variables = selectMatch[0]
            .replace(/SELECT\s+/i, '')
            .trim()
            .split(/\s+/)
            .map(v => v.replace('?', ''));
        
        return variables;
    } catch (e) {
        console.log('Could not extract variable order from query:', e);
        return null;
    }
}

function generateColumns(vars) {
    return vars.map(varName => ({
        title: varName.toUpperCase(),
        field: varName,
        sorter: "string",
        headerFilter: true
    }));
}

function transformData(bindings, vars) {
    return bindings.map(binding => {
        let row = {};
        vars.forEach(varName => {
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

    // Add chart cleanup function
    function cleanupChart() {
        const existingChart = Chart.getChart("myChart");
        if (existingChart) {
            existingChart.destroy();
        }
    }

    xhr.onload = function() {
        if (xhr.status >= 200 && xhr.status < 300) {
            try {
                var response = JSON.parse(xhr.responseText);
                console.log(response.message);
                
                // Check for null or empty response and throw error
                if (!response.message) {
                    cleanupChart(); // Cleanup chart before throwing error
                    throw new Error('Invalid query: No results returned. Please check your query syntax.');
                }
                
                // If there's empty array response (valid query but no matching results)
                if (response.message.length === 0) {
                    cleanupChart(); // Cleanup chart for empty results
                    table.setColumns([]);
                    table.setData([]);
                    document.getElementById('queryResults').innerHTML = 'Query executed successfully, but no results were found.';
                    return;
                }

                // Try to get order from query, fall back to response keys if needed
                let vars = extractVariableOrder(query);
                if (!vars || vars.length === 0) {
                    console.log('Falling back to response order');
                    vars = Object.keys(response.message[0]);
                }
                
                console.log('Variables order:', vars);

                var columns = generateColumns(vars);
                var data = transformData(response.message, vars);
                
                console.log('Transformed data:', data);
                console.log('Generated columns:', columns);
                
                table.setColumns(columns);
                table.setData(data);

                document.getElementById('downloadBtn').classList.replace("d-none", "d-md-flex");

                // Clear any previous error messages
                document.getElementById('queryResults').innerHTML = '';

                const queryType = determineQueryType(query, response.message);
                if (queryType) {
                    generateChart(response.message);
                } else {
                    // Clean up any existing chart if we don't have a visualization for this query
                    const existingChart = Chart.getChart("myChart");
                    if (existingChart) {
                        existingChart.destroy();
                    }
                }

            } catch (error) {
                console.error('Error processing SPARQL results:', error);
                // Don't show the Chart.js error to users, only show query-related errors
                const errorMessage = error.message.includes('Chart') 
                    ? 'Error processing query results.' 
                    : error.message;
                document.getElementById('queryResults').innerHTML = 'Please insert a valid query. ' + errorMessage;
                
                // Clear the table when there's an error
                table.setColumns([]);
                table.setData([]);
                
                // Hide download button if there's an error
                document.getElementById('downloadBtn').classList.replace("d-md-flex", "d-none");
                
                // Ensure chart is cleaned up on error
                cleanupChart();
            }
        } else {
            cleanupChart();
            document.getElementById('queryResults').innerHTML = 'Failed to load data: Status ' + xhr.status;
        }
    };

    xhr.onerror = function() {
        cleanupChart();
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
    query3: `SELECT ?humanReadableID ?admixtureName ?admixtureDensity
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
    query4: `SELECT DISTINCT ?humanReadableID ?compressiveStrength ?WaterCementRatio
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
    query5: `SELECT DISTINCT ?humanReadableID ?compressiveStrength ?elasticModulus
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
    query6: `SELECT DISTINCT ?humanReadableID ?E_Module ?WaterCementRatio
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

// map of query types to their visualization settings
const queryVisualizationSettings = {
    // Query with water-cement ratio (2 variables)
    'waterCementRatio': {
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
                responsive: true,
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const dataPoint = data[context.dataIndex];
                                return `${dataPoint.humanReadableID}: ${dataPoint.waterCementRatio}`;
                            }
                        }
                    }
                },
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Water-Cement Ratio'
                        }
                    }
                }
            }
        })
    },
    
    // Query with strength and modulus (3 variables)
    'strengthModulus': {
        type: 'scatter',
        config: (data) => ({
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Strength vs E-Module',
                    data: data.map(d => ({
                        x: parseFloat(d.compressiveStrength),
                        y: parseFloat(d.elasticModulus),
                        label: d.humanReadableID
                    })),
                    backgroundColor: 'rgba(75, 192, 192, 0.5)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    pointRadius: 6
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                const point = context.raw;
                                return `ID: ${point.label}, Strength: ${point.x}, E-Module: ${point.y}`;
                            }
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'Compressive Strength'
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'Elastic Modulus'
                        }
                    }
                }
            }
        })
    },
    
    // Query with admixture content (3 variables)
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

    'line': {
        type: 'scatter',  // Changed to scatter to allow custom x and y values
        config: (data) => ({
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'Compressive Strength vs W/C Ratio',
                    data: data.map(d => ({
                        x: parseFloat(d.compressiveStrength),
                        y: parseFloat(d.WaterCementRatio),
                        humanReadableID: d.humanReadableID  // Store ID for tooltip
                    })),
                    borderColor: 'rgb(75, 192, 192)',
                    backgroundColor: 'rgba(75, 192, 192, 0.5)',
                    pointRadius: 6,
                    showLine: true  // This will connect the points
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
                                    `WC Ratio: ${dataPoint.y}`,
                                    `Strength: ${dataPoint.x} MPa`
                                ];
                            }
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'Compressive Strength (MPa)'
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'Water-Cement Ratio'
                        }
                    }
                }
            }
        })
    },
    'line2': {
        type: 'scatter',  // Changed to scatter to allow custom x and y values
        config: (data) => ({
            type: 'scatter',
            data: {
                datasets: [{
                    label: 'E-Module vs W/C Ratio',
                    data: data.map(d => ({
                        x: parseFloat(d.E_Module),
                        y: parseFloat(d.WaterCementRatio),
                        humanReadableID: d.humanReadableID  // Store ID for tooltip
                    })),
                    borderColor: 'rgb(75, 192, 192)',
                    backgroundColor: 'rgba(75, 192, 192, 0.5)',
                    pointRadius: 6,
                    showLine: true  // This will connect the points
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
                                    `WC Ratio: ${dataPoint.y}`,
                                    `Strength: ${dataPoint.x}`
                                ];
                            }
                        }
                    }
                },
                scales: {
                    x: {
                        title: {
                            display: true,
                            text: 'E-Module'
                        }
                    },
                    y: {
                        title: {
                            display: true,
                            text: 'Water-Cement Ratio'
                        }
                    }
                }
            }
        })
    }
};

// Function to determine which visualization to use based on the query
function determineQueryType(query, data) {
    if (!data || data.length === 0) return null;
    
    const variables = Object.keys(data[0]);
    console.log(11111, variables)
    
    // Check for specific variable combinations
    if (variables.includes('waterCementRatio') && variables.length == 2) {
        console.log(11)
        return 'waterCementRatio';
    } else if(variables.includes('WaterCementRatio') && variables.includes('compressiveStrength')){
        console.log(22)
        return 'line';
    }else if(variables.includes('WaterCementRatio') && variables.includes('E_Module')){
        console.log(22)
        return 'line2';
    }else if (variables.includes('compressiveStrength') && variables.includes('elasticModulus')) {
        console.log(33)
        return 'strengthModulus';
    } else if (variables.includes('admixtureName')) {
        console.log(44)
        return 'admixture';
    } else if (variables.includes('cementType')) {
        console.log(55)
        return 'cement';
    }
    console.log(66)
    return null; // No visualization defined for this query type
}

function generateChart(data) {
    // Clean up existing chart
    const existingChart = Chart.getChart("myChart");
    if (existingChart) {
        existingChart.destroy();
    }
    
    // Get canvas element
    const ctx = document.getElementById('myChart');
    if (!ctx) return;
    
    // Determine query type and get visualization settings
    const queryType = determineQueryType(null, data);
    if (!queryType) {
        console.log('No visualization defined for this query type');
        return; // Don't create a chart if we don't have a visualization defined
    }
    
    const settings = queryVisualizationSettings[queryType];
    const chartConfig = settings.config(data);
    
    // Create new chart
    new Chart(ctx, chartConfig);
}