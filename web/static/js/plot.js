fetch('getJson')
    .then(response => response.json())
    .then(data => {
        // Extract the variables from the response
        const jsonList = data.json;

        // Log the variables
        console.log(jsonList);
    })
    .catch(error => {
        // Handle errors
        console.error('Fehler:', error);
    });

var trace1 = {
    x: [1, 2, 3, 4],
    y: [10, 15, 13, 17],
    mode: 'markers',
    type: 'scatter'
};

var trace2 = {
    x: [2, 3, 4, 5],
    y: [16, 5, 11, 9],
    mode: 'lines',
    type: 'scatter'
};

var data = [trace1, trace2];

Plotly.newPlot('myDiv', data);