function getMixtures() {
    // Get the mixtureId from the URL if it exists
    const urlParams = new URLSearchParams(window.location.search);
    const urlMixtureId = urlParams.get('mixtureId');

    const url = window.appConfig.urls.getMixtures;
    // Fetch the list of mixtures from the server
    fetch(url, {
        method: 'GET',
        headers: {
            'Content-Type': 'application/json',
        },
    })
    .then(response => {
        if (!response.ok) {
            throw new Error('Network response was not ok');
        }
        return response.json(); // Expect the response from the backend as JSON
    })
    .then(data => {
        // Process the response from the backend
        console.log(data); // Log the data to check its structure
        var select = document.getElementById('mixtures');
        // Loop over the mixtures and add each one as an option to the select
        for (var i = 0; i < data.mixtures.length; i++) {
            var option = document.createElement('option');
            option.value = data.mixtures[i].id;
            option.text = data.mixtures[i].name;
            select.appendChild(option);
            // If this mixture matches the one in the URL, select it
            if (data.mixtures[i].id === urlMixtureId) {
                option.selected = true;
                mixtureID = urlMixtureId;
                mixtureName = data.mixtures[i].name;
            }
        }
        // This line updates the select with Select2
        $(select).select2();
    })
    .catch(error => {
        console.error('There was a problem with your fetch operation:', error);
    });
}

$(document).ready(function() {
    getMixtures();
    $('#mixtures').select2();
});