function toggleSpecimenDimensions() {
    const shape = document.getElementById('SpecimenShape').value;
    const heightDiv = document.getElementById('SpecimenHeightDiv');
    const diameterDiv = document.getElementById('SpecimenDiameterDiv');
    const lengthDiv = document.getElementById('SpecimenLengthDiv');
    const widthDiv = document.getElementById('SpecimenWidthDiv');

    if (shape === 'cylindrical') {
        heightDiv.style.display = 'block';
        diameterDiv.style.display = 'block';
        lengthDiv.style.display = 'none';
        widthDiv.style.display = 'none';
    } else if (shape === 'cube' || shape === 'prism') {
        heightDiv.style.display = 'block';
        diameterDiv.style.display = 'none';
        lengthDiv.style.display = 'block';
        widthDiv.style.display = 'block';
    } else {
        heightDiv.style.display = 'none';
        diameterDiv.style.display = 'none';
        lengthDiv.style.display = 'none';
        widthDiv.style.display = 'none';
    }
}

