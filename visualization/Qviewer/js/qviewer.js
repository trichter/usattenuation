/* Qviewer | Copyright 2017 Tom Eulenfeld | Licensed MIT */
$(window).load(function(){

// http://stackoverflow.com/a/5077091
String.prototype.format = function () {
  var args = arguments;
  return this.replace(/\{(\d+)\}/g, function (m, n) { return args[n]; });
};

// http://stackoverflow.com/a/24094619
Array.prototype.SumArray = function (arr) {
    var sum = [];
    if (arr !== null && this.length == arr.length) {
        for (var i = 0; i < arr.length; i++) {
            sum.push(this[i] + arr[i]);
        }
    }
    return sum;
};
Array.prototype.DifArray = function (arr) {
    var sum = [];
    if (arr !== null && this.length == arr.length) {
        for (var i = 0; i < arr.length; i++) {
            sum.push(this[i] - arr[i]);
        }
    }
    return sum;
};

//http://stackoverflow.com/a/10284006
function zip(arrays) {
    return arrays[0].map(function(_,i){
        return arrays.map(function(array) {return array[i];});
    });
}

function log10(val) {
  return Math.log(val) / Math.LN10;
}

function pow10(val) {
  return Math.pow(10, val);
}

var freqs = ["1.5Hz", "2.1Hz", "3.0Hz", "4.2Hz", "6.0Hz", "8.5Hz", "12Hz", "17Hz", "slope"];
var freqsim = ["_01.5Hz", "_02.1Hz", "_03.0Hz", "_04.2Hz", "_06.0Hz", "_08.5Hz", "_12.0Hz", "_17.0Hz", "_slope"];
var image = "data/{0}{1}{2}{3}.png";

var station_size_ratio = 200;
var sta_index = -1;
function quan() {return $('#radioform input:radio:checked').val();}
function freqi() {return $("#slider").slider("value");}
function interp() {return $('#check3').is(':checked');}
function error() {return $('#check4').is(':checked');}
function sta() {return $('#check1').is(':checked');}
function noslope() {return (quan() == 'R') || (quan() == 'nobs') || (quan() == 'Qtot') || (quan() == 'B0');}
function noerror() {return (quan() == 'nobs') || (quan() == 'Qtot') || (quan() == 'B0');}
function derivative() {return (quan() == 'Qtot') || (quan() == 'B0');}

function update_slidertext(i) {
    if (typeof i === 'undefined'){i = freqi();}
    $(".slidertext").html(freqs[i]);
}
function update_image(i) {
    if (typeof i === 'undefined'){i = freqi();}
    if (quan() == "Nothing") {
        var src = image.format("us", "", "", "");
    } else {
        var s1 = "";
        var s2 = "";
        if (error()){s1 = "_error";}
        if (interp()){s2 = "_int";}
        var src = image.format(quan(), s1, freqsim[i], s2);
    }
    $("#image").attr("src", src);
}
function update_stations() {
    if (sta()) {
        plot1.setData([plotstations]);
    } else {
        plot1.unhighlight();
        plot1.setData([]);
    }
    plot1.draw();
}
function plot_stuff() {
    var sta = "";
    if ((sta_index > -1)) {
        sta = stations[0][sta_index]
    }    
    $(".labeltext").html(sta);
    if ((quan() == "Nothing") || (sta_index < 0) ||
            (derivative()) && (!results["Qsc"].hasOwnProperty(sta)) ||
            (!derivative()) && (!results[quan()].hasOwnProperty(sta))) {
        $(".xlabeltext").html("");
        $(".ylabeltext").html("");
        options2.grid.show = false;
        plot2 = $.plot(placeholder2, [], options2);
        return;
    }
    options2.grid.show = true;
    if (quan() == 'nobs') {
        var d = results[quan()][sta];
        options2.yaxis.min = null;
        options2.yaxis.max = null;        
        plotQ.points.errorbars = null;
        plotQ.data = zip([results.freqs.map(log10), d]);
        $(".ylabeltext").html("nobs");
    }
    else
    {
        if (derivative()) {
            var d1 = results["Qsc"][sta];
            var d2 = results["Qi"][sta];
            var d3 = d1.mean.map(pow10).SumArray(d2.mean.map(pow10)).map(log10);
            if (quan() == "Qtot") {
                var d = d3;
                var ymax = Math.max.apply(Math, d);
                var ymin = Math.min.apply(Math, d);
                var t = "- log10 Qtot";
            }
            else
            {
                var d = d1.mean.DifArray(d3).map(pow10);
                var ymax = 1;
                var ymin = 0;
                var t = "B0";
            }
            plotQ.points.errorbars = null;
            options2.yaxis.max = Math.ceil(2*ymax)/2;
            options2.yaxis.min = Math.floor(2*ymin)/2;           
            plotQ.data = zip([results.freqs.map(log10), d]);            
        }
        else
        {
            var d = results[quan()][sta];
            // setting null of errors to 0
            var i, n = d.error.length;
            for (i = 0; i < n; ++i) {
                if (d.error[i] == null) {
                    d.error[i] = 0;
                }
            }        
            var ymax = Math.max.apply(Math, d.mean.SumArray(d.error));
            options2.yaxis.max = Math.ceil(2*ymax)/2;
            var ymin = Math.min.apply(Math, d.mean.DifArray(d.error));
            options2.yaxis.min = Math.floor(2*ymin)/2;           
            plotQ.points.errorbars = "y";
            plotQ.data = zip([results.freqs.map(log10), d.mean, d.error]);
            var t = "log10 " + quan();
            if (!noslope()) {
                t = "- " + t;
            }
        }
        $(".ylabeltext").html(t);        
    }
    if (noslope()) {
        plot2 = $.plot(placeholder2, [plotQ], options2);
    } else {
        plotslope.data = [[0, d.intercept], [1.4, d.intercept + 1.4*d.slope]];
        plot2 = $.plot(placeholder2, [plotslope, plotQ], options2);
    }
    $(".xlabeltext").html("frequency (Hz)");
}

// Define slider and what happens when sliding
$("#slider").slider({
    min: 0,
    max: freqs.length - 1,    
    value: 0,
    slide: function(event, ui) {
        i = ui.value;
        if ((i == freqs.length - 1) && (noslope() || error())){
            return false;
        }
        if (!noslope()) {      
            $("#check4").prop("disabled", i == freqs.length - 1);
            $("#check4").button("refresh");
        }
        update_slidertext(i);
        update_image(i);
    },
    change: function(event, ui) {
        update_slidertext();
        update_image();
    }
});

// Radiobutton action
$("#radio").buttonset();
$("#radioform input:radio").change( function(){
    if (noerror()) {
        $("#check4").prop("checked", false);
    }
    $("#check4").prop("disabled", (noerror()) || ((freqi() == freqs.length - 1) && (!noslope())));
    $("#check4").button("refresh");
    if ((freqi() == freqs.length - 1) && (noslope() || error())){
        $("#slider").slider('value', freqi() - 1);
        update_slidertext();
    }
    update_image();
    plot_stuff();
});

// Checkbutton action
$("#check").buttonset();
$('#checkform input:checkbox').change(function () {
    switch ($(this).val()) {
        case "interpolate":
            update_image();
            break;
        case "error":           
            update_image();
            break;
        case "stations":
            update_stations();
            break;
    }
});

// Define plots
var placeholder1 = $("#placeholder1");
var placeholder2 = $("#placeholder2");

var plotstations = {
    data: stations[1],
    color: "#000000",
    points: {
        fill: false, lineWidth: 1,
        radius: placeholder1.width() / station_size_ratio
    },
    clickable: false,
    hoverable: true,
    highlightColor: "#dba255",    
};
var plotQ = {
    data: null,
    color: "#000000",
    lines: {show: false},
    points: {
        show: true, fill: false, lineWidth: 1, errorbars: "y",
        yerr: {show: true, upperCap: "-", lowerCap: "-", color: "#808080"}
    }
};
var plotslope = {
    data: null,
    color: "#808080",
    points: {show: false},
    lines: {show: true, zero:false}
};   
var options1 = {
    grid: {       
        show: false,
        margin: 0, 
        clickable: true,
        hoverable: true     
    },
    series: {
        lines: {show: false},
        points: {show: true}
    },
    xaxis: {min: 0, max: 1},
    yaxis: {min: 0, max: 1}     
};
var options2 = {
    grid: {show: false},
    yaxis: {},
    xaxis: {
        min: 0, max: 1.406,
        ticks: [[0, "1"], [0.301, ""], [0.477, ""], [0.602, ""], [0.699, ""],
                [0.778, ""], [0.845, ""], [0.903, ""], [0.954, ""], [1, "10"],
                [1.30, "20"]]
    }
};

var plot1 = $.plot(placeholder1, [], options1);
var plot2 = $.plot(placeholder2, [], options2);

// Define station tooltip
$("<div id='stationtooltip'></div>").css({
    position: "absolute",
    display: "none",
    border: "1px solid #fdd",
    "border-radius": "5px",
    padding: "2px",
    "background-color": "#fee",
    "z-index": 3,
    opacity: 0.80
}).appendTo("body");

// Display tooltip if stations are hovered
// Plot station results if station is clicked in the map
// Hide station if Ctr clicked
// http://stackoverflow.com/questions/7645830/flot-detect-left-click-and-right-click
var currentItem = null;
placeholder1.bind("plothover", function (event, pos, item) {
    if (item) {
        currentItem = item;
        $("#stationtooltip").html(stations[0][item.dataIndex])
        .css({top: item.pageY+10, left: item.pageX+10})
        .fadeIn(200);
    } else {
        currentItem = null;
        $("#stationtooltip").hide();
    }
});

placeholder1.mousedown(function(event) {
    var item = currentItem;
    if (item) {            
        if (event.ctrlKey) {
            if (item.dataIndex == sta_index) {
                plot1.unhighlight();
                sta_index = -1;
            }
            stations[1][item.dataIndex] = [-1, -1];            
            plot1.setData([plotstations]);
            plot1.draw();            
        } else {
            plot1.unhighlight();
            plot1.highlight(item.series, item.datapoint);
            sta_index = item.dataIndex;
        }
        plot_stuff();        
    }
});

/*placeholder1.bind("plotclick", function (event, pos, item) {
    if (item) {        
        plot1.unhighlight();
        plot1.highlight(item.series, item.datapoint);
        sta_index = item.dataIndex;
        plot_stuff();        
    }
});*/

// Adapt circle size of stations when resizing
placeholder1.resize(function () {
    plotstations.points.radius = $(this).width() / station_size_ratio;
    if (sta()) {plot1.setData([plotstations]);}
});

// Define draggables
$("#drag1").resizable({
    aspectRatio: true,
    maxWidth: 1500,
    minWidth: 100
});
$("#drag1").draggable();
$("#drag2").resizable({
    aspectRatio: false,
    maxWidth: 500,
    minWidth: 50,
    maxHeight: 500,
    minHeight: 50
});
$("#drag2").draggable();
//$("#controls").draggable({cancel: ".radio"});

// Use qtip2 for nice tooltips
$('[title]').qtip({
    show: {delay: 1000},
    style: {tip: false, classes: "qtip-rounded"},
    position: {at: 'bottom center', adjust: {x: 10, y: 5}}
});

// Load first image, slider text, plot stations
update_slidertext();
update_image();
update_stations();

});
