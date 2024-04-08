document.addEventListener('DOMContentLoaded', function() {
    // Funktion zum Ausführen der SPARQL-Abfrage
    function executeSparqlQuery(query) {
        // Erstellen eines neuen XMLHttpRequest
        var xhr = new XMLHttpRequest();
        // Konfigurieren der Anfrage
        xhr.open('POST', '/query', true); // oder GET, abhängig von Ihrem Endpunkt
        xhr.setRequestHeader('Content-Type', 'application/x-www-form-urlencoded'); // Für POST-Anfragen

        // Setzen eines Handlers für die Antwort
        xhr.onload = function() {
            if (xhr.status >= 200 && xhr.status < 300) {
                // Erfolg: Ergebnisse im angegebenen Container anzeigen
                document.getElementById('queryResults').innerHTML = xhr.responseText;
            } else {
                // Fehler: Fehlermeldung im Container anzeigen
                document.getElementById('queryResults').innerHTML = 'Fehler bei der Ausführung der Abfrage.';
            }
        };

        // Senden der Anfrage mit Daten
        xhr.send('query=' + encodeURIComponent(query)); // Stellen Sie sicher, dass dies Ihrem Endpunkt-Format entspricht
    }

    // Formular-Submit-Event-Handler
    document.getElementById('sparqlForm').addEventListener('submit', function(e) {
        e.preventDefault(); // Standard-Formular-Submit unterbrechen
        var query = document.getElementById('sparqlQuery').value;

        // Führen Sie die Abfrage mit der definierten Funktion aus
        executeSparqlQuery(query);
    });
    // Event-Handler für den Button, um den vorgefertigten Text einzufügen
    document.getElementById('insertTextBtn').addEventListener('click', function() {
        var predefinedText = 'SELECT ?s ?p ?o WHERE { ?s ?p ?o }';
        document.getElementById('sparqlQuery').value = predefinedText;
    });
});

//js for upload data page
// Globale Variable
var mixtureID = null;
var institute = 'BAM';

function generateUUIDv4() {
return ([1e7]+-1e3+-4e3+-8e3+-1e11).replace(/[018]/g, c =>
    (c ^ crypto.getRandomValues(new Uint8Array(1))[0] & 15 >> c / 4).toString(16)
);
}

function uploadData(type, buttonID) {
    if (type === 'Mixture') {
    mixtureID = generateUUIDv4(); // Generiere eine neue UUID v4 und weise sie `mixtureID` zu
        console.log(`Neue Mixture ID: ${mixtureID}`); // Zeige die generierte ID an
    } else {
    console.log(`Type ist nicht Mixture, keine ID generiert. Mixture ID ist: ${mixtureID}`);
    }
    var fileInput = document.getElementById(buttonID);
    // Check the file format (extension)
    const allowedFormats = ['xlsx', 'xls', 'csv', 'dat', 'txt'];
    const fileExtension = fileInput.files[0].name.split('.').pop().toLowerCase();

    console.log(fileExtension)
    if (!allowedFormats.includes(fileExtension)) {
        alert('Invalid file. Please select a xlsx, xls, csv, dat or txt file.');
        return;
    }
    
    var formData = new FormData();
    formData.append('file', fileInput.files[0]); // Fügt die Datei hinzu
    formData.append('type', type); // Fügt den übergebenen Typ hinzu
    formData.append('Mixture_ID', mixtureID); // Übergibt die Mixture ID

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
        alert(data.message); // Zeigt eine Nachricht vom Server an
    })
    .catch((error) => {
        console.error('Fehler beim Hochladen:', error);
    });
}

function submitMixture() {
    var textInput = document.getElementById('textInput').value;
    var apiUrl = '/search'; // Pfad zur Flask-Route, relativ zur Basis-URL der Website

    // Macht den Aufruf an das Backend
    fetch(apiUrl, {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
        },
        body: JSON.stringify({ mixtureName: textInput }), // Sendet den Namen der Mischung als JSON
    })
        .then(response => {
        if (!response.ok) {
            throw new Error('Network response was not ok');
        }
        return response.json(); // Erwartet die Antwort des Backends als JSON
    })
        .then(data => {
        // Verarbeitet die Antwort des Backends
        var responseDiv = document.getElementById('mixresponseText');
        responseDiv.style.display = 'block'; // Macht den Antwort-Text sichtbar
        responseDiv.classList = 'text-danger';
        responseDiv.innerHTML = data.message; // Setzt die Antwort des Backends in das Div
        
        // Speichert die mixtureID als Variable, ohne sie anzuzeigen
        if (data.mixtureID) {
            mixtureID = data.mixtureID; // Speichert den Wert in der Variable
            responseDiv.classList = 'text-success';
            console.log('Gespeicherte mixtureID:', mixtureID); // Optional: Zur Überprüfung in der Konsole ausgeben
}
    })
        .catch(error => {
        console.error('There was a problem with your fetch operation:', error);
        var responseDiv = document.getElementById('mixresponseText');
        responseDiv.style.display = 'block';
        responseDiv.classList = 'text-danger';
        responseDiv.innerHTML = 'Fehler beim Abrufen der Daten'; // Zeigt eine Fehlermeldung an
    });
}

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