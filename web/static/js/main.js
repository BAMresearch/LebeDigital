//js for upload data page
// Globale Variable
var mixtureID = null;
var mixtureName = null;
var institute = 'BAM';

function showMixtureName() {
    if (mixtureID !== null) {
        var elements = document.getElementsByClassName('show-mixture-id');
        for (var i = 0; i < elements.length; i++) {
            elements[i].innerHTML = "Mixture: "+ mixtureName;
        }
    }    
}

function hideMixtureId() {
    if (mixtureID == null) {
        var elements = document.getElementsByClassName('show-mixture-id');
        for (var i = 0; i < elements.length; i++) {
            elements[i].innerHTML = "";
        }
    }    
}

function generateUUIDv4() {
return ([1e7]+-1e3+-4e3+-8e3+-1e11).replace(/[018]/g, c =>
    (c ^ crypto.getRandomValues(new Uint8Array(1))[0] & 15 >> c / 4).toString(16)
);
}

function uploadData(type, fileID, urlID, label) {
    if (type === 'Mixture') {
        mixtureID = generateUUIDv4(); // Generiere eine neue UUID v4 und weise sie `mixtureID` zu
        console.log(`Neue Mixture ID: ${mixtureID}`); // Zeige die generierte ID an
    } else {
        console.log(`Type ist nicht Mixture, keine ID generiert. Mixture ID ist: ${mixtureID}`);
    }
    
    var formData = new FormData();
    formData.append('type', type); // Fügt den übergebenen Typ hinzu
    formData.append('Mixture_ID', mixtureID); // Übergibt die Mixture ID

     // Check if the input is a file or a URL
    var fileLabel = document.getElementById(label).textContent;
    if (fileLabel.startsWith('http')) {
        // It's a URL, append it to the form data
        var urlInput = document.getElementById(urlID).value;
        formData.append('url', urlInput);

        // Extract the filename from the URL
        var url = new URL(urlInput);
        var pathname = url.pathname;
        var filename = pathname.substring(pathname.lastIndexOf('/') + 1);
        mixtureName = filename;

    } else {
        // It's a file
        var fileInput = document.getElementById(fileID);
        for (var i = 0; i < fileInput.files.length; i++) {
           formData.append('file' + i, fileInput.files[i]); // Fügt jede Datei hinzu
        }
        //formData.append('file', fileInput.files[0]); // Fügt die Datei hinzu
        mixtureName = fileLabel;
    }
    showMixtureName();

    // Log each key-value pair in formData
    for (var pair of formData.entries()) {
        console.log(pair[0] + ', ' + pair[1]);
    }

    fetch('/dataUpload', {
        method: 'POST',
        body: formData, // Sendet FormData, das Datei und Typ enthält
    })
    .then(response => {
        if (!response.ok) {
            throw new Error('Netzwerkantwort war nicht ok');
        }
        return response.json();
    })
    .then(data => {
        if (data.status === 200) {
            document.getElementById("message").innerHTML = data.message
            const toastLiveExample = document.getElementById('liveToast')
            const toastBootstrap = bootstrap.Toast.getOrCreateInstance(toastLiveExample)
            toastBootstrap.show()  
        } else {
            document.getElementById("error-message").innerHTML = data.message
            const toastLiveExample = document.getElementById('liveToastError')
            const toastBootstrap = bootstrap.Toast.getOrCreateInstance(toastLiveExample)
            toastBootstrap.show()
        }
    })
    .catch((error) => {
        console.error('Fehler beim Hochladen:', error);
    });
}

// This function is called when the upload button is clicked
function clearFileInput(fileID) {
    document.getElementById(fileID).value = '';
}

// This function is called when multiple files are selected
function onFileSelected(event, fileLabel) {
    var files = event.target.files; // Get the list of files
    var fileNames = []; // Initialize an array to store the names of valid files

    // Check the file format (extension) for each file
    const allowedFormats = ['xlsx', 'xls', 'csv', 'dat', 'txt', 'json', 'xml'];
    for (var i = 0; i < files.length; i++) {
        var fileName = files[i].name;
        const fileExtension = fileName.split('.').pop().toLowerCase();

        if (!allowedFormats.includes(fileExtension)) {
            alert('Invalid file: ' + fileName + '. Please select a xlsx, xls, csv, dat, txt, json, or xml file.');
            continue; // Skip this file and move to the next one
        }
        fileNames.push(fileName); // Add the valid file name to the array
    }

    // Display the names of valid files
    document.getElementById(fileLabel).textContent = fileNames.join(', ');
    //event.target.value = '';  // Reset the value
}


// This function is called when the "ok" button in the modal is clicked
function onUrlEntered(urlInput,fileLabel) {
    var url = document.getElementById(urlInput).value;
    if (isValidURL(url)) {
        document.getElementById(fileLabel).textContent = url;
    } else {
        alert("Please enter a valid URL.");
    }
    var modalElement = document.getElementById('exampleModal');
    var modalInstance = bootstrap.Modal.getInstance(modalElement);
    modalInstance.hide();
}

// This function checks if a URL is valid
function isValidURL(string) {
    var res = string.match(/(http|https):\/\/(\w+:{0,1}\w*@)?(\S+)(:[0-9]+)?(\/|\/([\w#!:.?+=&%@!\-\/]))?/);
    return (res !== null)
};


$(document).ready(function() {
    getMixtures();
    $('#mixtures').select2();
});

function getMixtures() {
    // Fetch the list of mixtures from the server
    fetch('/get-mixtures', {
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
        }
        // This line updates the select with Select2
        $(select).select2();
    })
    .catch(error => {
        console.error('There was a problem with your fetch operation:', error);
    });
}

// Update the mixtureId variable and the label when a mixture is selected
$('#mixtures').on('change', function() {
    mixtureID = this.value;
    mixtureName = $(this).find('option:selected').text();
    showMixtureName();
});


// Nächste Seite
function nextSiteOne(page) {
    // Dictionary für seiten Namen
    let pageMap = new Map();

    pageMap.set(1, 'mainContent');
    pageMap.set(2, 'newContent');
    pageMap.set(3, 'last');
    
    // Verbirgt den aktuellen Hauptinhalt
    document.getElementById(pageMap.get(page)).style.display = 'none';

    // Zeigt den neuen Inhalt an
    document.getElementById(pageMap.get(page + 1)).style.display = 'block';
}

// previous page
function previousSiteOne(page) {
    // Dictionary für seiten Namen
    let pageMap = new Map();

    pageMap.set(1, 'mainContent');
    pageMap.set(2, 'newContent');
    pageMap.set(3, 'last');

    document.getElementById(pageMap.get(page)).style.display = 'block';

    document.getElementById(pageMap.get(page + 1)).style.display = 'none';
}

function toggleSections() {
    const existingMixture = document.getElementById("existingMixture");
    const mixtureDetails = document.getElementById("mixtureDetails");
    const fileUpload = document.getElementById("fileUpload");

    if (existingMixture.checked) {
        mixtureDetails.style.display = "block";
        fileUpload.style.display = "none";
    } else {
        mixtureDetails.style.display = "none";
        fileUpload.style.display = "block";
    }
}