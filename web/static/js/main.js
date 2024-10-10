//js for upload data page
// Global Variable
var mixtureID = null;
var mixtureName = null;
var institute = 'BAM';

function enableUploadButton() {
    document.querySelectorAll('.uploadbtn').forEach(button => {
        button.disabled = false;
    });
}

function disableUploadButton() {
    document.querySelectorAll('.uploadbtn').forEach(button => {
        button.disabled = true;
    });
}


// Function to display the selected mixture name on the UI
function showMixtureName() {
    console.log(mixtureName)
    if (mixtureID !== null) {
        var elements = document.getElementsByClassName('show-mixture-id');
        for (var i = 0; i < elements.length; i++) {
            elements[i].innerHTML = "Mixture: "+ mixtureName;
        }
    }    
}

// Function to remove the mixture ID and clear related UI elements
function removeMixtureId() {
    mixtureID = null
    mixtureName = null
    var elements = document.getElementsByClassName('show-mixture-id');
    for (var i = 0; i < elements.length; i++) {
        elements[i].innerHTML = "";
    }   
}


// Function to generate unique id
function generateUUIDv4() {
return ([1e7]+-1e3+-4e3+-8e3+-1e11).replace(/[018]/g, c =>
    (c ^ crypto.getRandomValues(new Uint8Array(1))[0] & 15 >> c / 4).toString(16)
);
}


function uploadData(type, fileID, urlID, label) {
    if (type === 'Mixture') {
        mixtureID = generateUUIDv4(); 
        console.log(`New Mixture ID: ${mixtureID}`); 
    } else {
        console.log(`Type ist not Mixture, No ID generated. Mixture ID is: ${mixtureID}`);
        if (mixtureID == null) {
            // Show error if no mixture is selected
            $('#mixtureWarningModal').modal('show');
            return;
        }
    }

    var formData = new FormData();
    formData.append('type', type);
    formData.append('Mixture_ID', mixtureID);

    // Check if the input is a file or a URL
    var fileLabel = document.getElementById(label).textContent;
    if (fileLabel.startsWith('http')) {
        // Handle URL input
        var urlInput = document.getElementById(urlID).value;
        formData.append('url', urlInput);

        // Extract the filename from the URL
        var url = new URL(urlInput);
        var pathname = url.pathname;
        var filename = pathname.substring(pathname.lastIndexOf('/') + 1);
        if (type === 'Mixture') {
            mixtureName = filename;
        }

    } else {
        // Handle file input
        var fileInput = document.getElementById(fileID);
        for (var i = 0; i < fileInput.files.length; i++) {
           formData.append('file' + i, fileInput.files[i]); 
        }
        if (type === 'Mixture') {
            mixtureName = fileInput.files[0].name; 
        }
    }
    showMixtureName();

    fetch('/dataUpload', {
        method: 'POST',
        body: formData,
    })
    .then(response => {
        if (!response.ok) {
            throw new Error('Network was not Okay!');
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
            if(type === 'Mixture' && data.status === 409){
                removeMixtureId()
            }
            document.getElementById("error-message").innerHTML = data.message
            const toastLiveExample = document.getElementById('liveToastError')
            const toastBootstrap = bootstrap.Toast.getOrCreateInstance(toastLiveExample)
            toastBootstrap.show()
        }
        disableUploadButton();
    })
    .catch((error) => {
        console.error('Failed to upload:', error);
        document.getElementById("error-message").innerHTML = "Failed to upload!"
        const toastLiveExample = document.getElementById('liveToastError')
        const toastBootstrap = bootstrap.Toast.getOrCreateInstance(toastLiveExample)
        toastBootstrap.show()
    });
}

// Function to restrict user clicking upload btn more than once
function clearFileInput(fileID) {
    document.getElementById(fileID).value = '';
    disableUploadButton();
}

// This function is called when multiple files are selected
function onFileSelected(event, fileLabel) {
    var files = event.target.files; // Get the list of files
    var fileNames = []; // Initialize an array to store the names of valid files

    // Check the file format (extension) for each file
    var allowedFormats = ['json']
    if(fileLabel == 'fileLabel1'){
        // for mixture
        allowedFormats.push('xls', 'xlsx');
    }
    else if(fileLabel == 'fileLabel2'){
        // for comSt
        allowedFormats.push('dat')
    }
    else{
        allowedFormats.push('xml')
    }
    
    for (var i = 0; i < files.length; i++) {
        var fileName = files[i].name;
        const fileExtension = fileName.split('.').pop().toLowerCase();

        if (!allowedFormats.includes(fileExtension)) {
            alert('Sorry! The file type is not supported!');
            continue; // Skip this file and move to the next one
        }
        fileNames.push(fileName); // Add the valid file name to the array
    }

    // Display the names of valid files
    document.getElementById(fileLabel).textContent = fileNames.join(', ');
    enableUploadButton();
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

    // Hide the url input modal
    var modalId = ['exampleModal', 'exampleModal2', 'exampleModal3'];
    modalId.forEach(function(modalId) {
        var modalElement = document.getElementById(modalId);
        var modalInstance = bootstrap.Modal.getInstance(modalElement);
        if (modalInstance) {
            modalInstance.hide();
        }
    });
    enableUploadButton();
}

// This function checks if a URL is valid
function isValidURL(string) {
    var res = string.match(/(http|https):\/\/(\w+:{0,1}\w*@)?(\S+)(:[0-9]+)?(\/|\/([\w#!:.?+=&%@!\-\/]))?/);
    return (res !== null)
};


// Update the mixtureId variable and the label when a mixture is selected
$('#mixtures').on('change', function() {
    mixtureID = this.value;
    mixtureName = $(this).find('option:selected').text();
    showMixtureName();
});


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

// Function to redirect user to the mixture upload section
function goToMixtureUpload() {
    $('#mixtureWarningModal').modal('hide');
    // Scroll to and open the mixture accordion
    $('#accordionFlushExample .accordion-item:first-child .accordion-button').click();
    $('html, body').animate({
        scrollTop: $("#accordionFlushExample").offset().top
    }, 1000);
}

// Redirect user to the mixture form
function GoToMixtureForm() {
    // Unselect all radio buttons
    const radioButtons = document.querySelectorAll('input[name="mixtureOption"]');
    radioButtons.forEach(radio => radio.checked = false);

    // Redirect
    window.location.href = "{{ url_for('new_mixture') }}";  
}

// Function to check if a mixture is selected before redirecting to another page
function checkMixtureAndRedirect(targetPage) {
    if (mixtureID == null) {
        $('#mixtureWarningModal').modal('show');
    } else {
        redirectToPage(targetPage);
    }
}


// Function to handle page redirection with mixture ID as a query parameter
function redirectToPage(page) {
    let url = new URL(page, window.location.origin);
    url.searchParams.append('mixtureId', mixtureID);
    window.location.href = url.toString();
}

// Redirect user to the compressive strength form
function GoToComStForm() {
    checkMixtureAndRedirect('/new_compressive_strength');
}

// Redirect user to the E-Module form
function GoToEModuleForm() {
    checkMixtureAndRedirect('/new_emodule');
}