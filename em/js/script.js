function getValues() {
    x = document.getElementById("rangevalue").value;
    v0js = document.getElementById("rangevalue2").value;
    d = document.getElementById("rangevalue3").value;
}

setInterval(getValues,100);