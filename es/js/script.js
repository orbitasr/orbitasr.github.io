// Função que é executada sempre que a página é aberta
function init() {
    document.getElementById('ico1').classList.add('select');
    document.getElementById('ico1b').classList.add('select');
    if (localStorage.getItem("mode") == "dark") {
        dark();
    }
    else {
        localStorage.setItem("mode", "light")
    }
}

function dark() {
    var element = document.body;
    element.classList.add("dark-mode");
}

function change() {
    var element = document.body;
    element.classList.toggle("dark-mode");
    if (localStorage.getItem("mode") == "dark") {
        localStorage.setItem("mode", "light");
    } else {
        localStorage.setItem("mode", "dark");
    }
}


function show(shown, hide) {
    document.getElementById(shown).style.display='block';
    document.getElementById(hide).style.display='none';
}

function openModal(a,b) {
    document.getElementById(a).style.display='block';
    document.getElementById(b).style.display='block';
}

function hide(a,b){
    document.getElementById(a).style.display='none';
    document.getElementById(b).style.display='none';
}

function nav(a,b,c,d,e,f){
    show(a,b);
    document.getElementById(c).style.color='#FFC700';
    document.getElementById(d).style.color='#000000';
    document.getElementById(e).classList.add('select');
    document.getElementById(f).classList.remove('select');  
}

function getValues() {
    m = document.getElementById("rangevalue").value;
    e = document.getElementById("number").value;
    o = document.getElementById("rangevalue2").value;
    mode = localStorage.getItem("mode");
    D1 = document.getElementById("rangevalue3").value;
    D2 = document.getElementById("number2").value;
}

setInterval(getValues,100);

function mostraPasso3() {
    var orbitas = pyscript.interpreter.globals.get('orbitas');
    if (orbitas) {
        document.getElementById("parte3").style.display='block';
    }
}

setInterval(mostraPasso3,100);

function limpaDiv(show){
    document.getElementById(show).style.display='block';
    titulo = document.getElementById("titulo");
    titulo.style.display="block";
    info = document.getElementById("valores");
    info.style.display='block';
    info.firstChild.remove();
    var graph = document.getElementById("graph");
    graph.firstChild.remove();
    graph.firstChild.remove();
    return false
}

function limpaDiv2() {
    var titulo2 = document.getElementById("titulo2");
    titulo2.style.display="block";
    var graph = document.getElementById("graph2");
    graph.firstChild.remove();
    var info = document.getElementById("infos");
    info.firstChild.remove();
    
}

function limpaDiv3() {
    var graph = document.getElementById("graph3");
    graph.firstChild.remove();
}
