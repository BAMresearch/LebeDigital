var jsonList = null;

fetch('getJson')
    .then(response => response.json())
    .then(data => {
        // Extract the variables from the response
        jsonList = data.json;

        // Log the variables
        console.log(jsonList);
    })
    .catch(error => {
        // Handle errors
        console.error('Fehler:', error);
    });


function create_plot() {
    x_name = document.getElementById("input1").value;
    x_type = document.getElementById("type1").value;
    y_type = document.getElementById("type2").value;
    y_name = document.getElementById("input2").value;

    console.log(x_type);
    x = [];
    y = [];

    if (x_type == "Mixture" && y_type == "Mixture") {
        for (var i = 0; i < jsonList[x_type].length; i++) {
        // every element of json List
            var obj = jsonList[x_type][i];
            if (x_name in obj) {
                x.push(obj[x_name]);
            }
            
            if (y_name in obj) {
                y.push(obj[y_name]);
            }
        }
    } else {
        for (var i = 0; i < jsonList[x_type].length; i++) {
            for (var j = 0; j < jsonList[x_type][i].length; j++) {
                console.log(jsonList[x_type][i][j]);
                var obj = jsonList[x_type][i][j];
                if (x_name in obj) {
                    x.push(obj[x_name]);
                }
                if (y_name in obj) {
                    y.push(obj[y_name]);
                }
            }
        }
    }

    console.log(x);
    console.log(y);

    var trace1 = {
        x: x,
        y: y,
        mode: 'markers',
        type: 'scatter'
    };

    var trace2 = {
        x: [2, 3, 4, 5],
        y: [16, 5, 11, 9],
        mode: 'lines',
        type: 'scatter'
    };

    var data = [trace1];

    Plotly.newPlot('myDiv', data);
}