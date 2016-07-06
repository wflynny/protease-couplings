var inited = false,
    width = 700,
    height = 700,
    innerRadius = Math.min(width, height) * 0.43,
    outerRadius = innerRadius * 1.08

//create number formatting functions
var formatPercent = d3.format("%");
var numberWithCommas = d3.format("0,f");

//create the arc path data generator for the groups
var arc = d3.svg.arc()
    .innerRadius(innerRadius)
    .outerRadius(outerRadius);

//create the chord path data generator for the chords
var path = d3.svg.chord()
    .radius(innerRadius);

var last_layout; //store layout between updates
var protein, wt;
var g, svg;
var chordOpac = 0.75,
    arcOpac = 0.95;
var red = "#ca0020",
    blue = "#0571b0",
    orange = "#f4a582",
    darkorange = "#d95f02";
var chordFill = d3.scale.ordinal().domain(d3.range(2)).range([blue, red]);

function init() {
    svg = d3.select("#chart_placeholder").append("svg")
            .attr("width", width)
            .attr("height", 0)
    svg.transition().duration(500)
            .attr("height", height)

    g = svg.append("g")
            .style("opacity", 1e-6)
            .attr("transform",
                  "translate(" + width / 2 + "," + height / 2 + ")scale(0.25)");
    g.transition().delay(500).duration(750)
            .attr("id", "circle")
            .style("opacity", 1)
            .attr("transform",
                  "translate(" + width / 2 + "," + height / 2 + ")scale(1)");
    g.append("circle")
        .attr("r", outerRadius*1.05);

    /*d3.select("#info").selectAll(".hide")
        .transition().duration(200)
          .style("opacity", 0)
        .remove();*/

    wt = "PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMNLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
    d3.csv("data/protease-protein.csv", function(error, proteinData) {
        if (error) {alert("Error reading file: ", error.statusText); return; }
        protein = proteinData;
    });
    setTimeout(function() { window.scrollBy(0, 105); }, 400);

}

function initButtons() {
  d3.select("#wildtype-button")
    .transition().duration(1500)
      .style("opacity", 1);
  d3.select("#mutant1-button")
    .transition().duration(1500)
      .style("opacity", 1);
  d3.select("#mutant2-button")
    .transition().duration(1500)
      .style("opacity", 1);
  d3.select("#mutant3-button")
    .transition().duration(1500)
      .style("opacity", 1);
}

function hideZeros(d) {
    if (d.source.index == d.target.index)
   { return 1e-6; } else { return chordOpac; }
}

function getDefaultLayout() {
    return d3.layout.chord()
    .padding(0.03)
    .sortSubgroups(d3.descending)
    .sortChords(d3.ascending);
}

function updateChords( datasetURL, seqURL ) {

    d3.json(datasetURL, function(error, matrix) {

    if (error) {alert("Error reading file: ", error.statusText); return; }

    var values = [],
        interactions = [];

    var sequence = matrix.sequence;
    matrix.matrix.forEach(function(row) {
      values.push(row.map(Math.abs));
      interactions.push(row.map(Math.sign));
    });
    values.forEach(function(vals, index, arr) {
      if (vals.every(function(v) { return v == 0; })) {
        values[index][index] = 30;
      }
    });

    /* Compute chord layout. */
    layout = getDefaultLayout(); //create a new layout object
    layout.matrix(values);

    /* Create/update "group" elements */
    var groupG = g.selectAll("g.group")
        .data(layout.groups(), function (d) {
            return d.index;
        });

    groupG.exit()
        .transition()
            .duration(1500)
            .style("opacity", 1e-6)
            .remove(); //remove after transitions are complete

    var newGroups = groupG.enter().append("g")
        .attr("class", "group");
    //the enter selection is stored in a variable so we can
    //enter the <path>, <text>, and <title> elements as well

    //Create the title tooltip for the new groups
    newGroups.append("title");

    //Update the (tooltip) title text based on the data
    groupG.select("title")
        .text(function(d, i) {
            return "Position "+protein[i].position+sequence[i];
        });

    //create the arc paths and set the constant attributes
    //(those based on the group index, not on the value)
    newGroups.append("path")
        .attr("id", function (d) {
            return "group" + d.index;
        })
        .style("fill", function (d) { return arcMutColor(sequence, d) });
            //return protein[d.index].color;
//        });

    //update the paths to match the layout
    groupG.select("path")
        .transition()
            .duration(250)
            .style("opacity", 0.5*arcOpac)
        .attrTween("d", arcTween( last_layout ))
            .transition().duration(250).style("opacity", arcOpac)
            .style("fill", function (d) { return arcMutColor(sequence, d) })
        ;

    function arcMutColor(s, d) {
      var ch = s[d.index];
      return wt[d.index] != ch ? orange : "#AAAAAA";
    }
    function textMutColor(s, d) {
      var ch = s[d.index];
      return wt[d.index] != ch ? darkorange : "black";
    }

    //create the group labels
    var text = newGroups.append("svg:text");
    text.attr("xlink:href", function (d) { return "#group" + d.index; })
        .attr("dy", ".35em");
    text.append("svg:tspan")
        .text(function (d) { return protein[d.index].position });
    text.append("svg:tspan")
        .attr("class", "aa")
        .style("font-weight", "bold")
        .style("fill", function (d) { return textMutColor(sequence, d) })
        .text(function (d) { return sequence[d.index] });

    //insure labels are correct
    g.selectAll("g.group text tspan.aa")
        .style("fill", function (d) { return textMutColor(sequence, d) })
        .text(function (d) { return sequence[d.index] });

    function labelTransform(d) {
      d.angle = (d.startAngle + d.endAngle) / 2;
      degAngle = d.angle * 180 / Math.PI;
      lt180 = (90 - degAngle)/2;
      gt180 = (-90 - degAngle)/2;
      x = degAngle > 180 ? outerRadius + 6 :
        outerRadius + 6*Math.cos(lt180*Math.PI/180);
      y = degAngle > 180 ? +6*Math.sin(gt180*Math.PI/180) :
        -6*Math.sin(lt180*Math.PI/180);

      return "rotate(" + (degAngle-90) + ")" +
          " translate(" + x + ',' + y + ")" +
          (degAngle > 180 ? " rotate("+gt180+")" :
           " rotate("+lt180+")");
    }
    function flipLabel(d) {
      return d.angle > Math.PI ? "end" : "begin";
    }

    //position group labels to match layout
    groupG.select("text")
        .transition().duration(500)
        .attr("transform", labelTransform)
        .attr("text-anchor", flipLabel);

    /* Create/update the chord paths */
    var chordPaths = g.selectAll("path.chord")
        .data(layout.chords(), chordKey );

    //create the new chord paths
    var newChords = chordPaths.enter()
        .append("path")
        .attr("class", "chord");

    // Add title tooltip for each new chord.
    newChords.append("title");

    // Update all chord title texts
    chordPaths.select("title")
        .text(function(d) {
            return ["Position ",
                    protein[d.source.index].position+sequence[d.source.index],
                    (interactions[d.source.index][d.target.index] < 0 ?
                     " positively " : " negatively "), "coupled with position ",
                    protein[d.target.index].position+sequence[d.target.index]
                    ].join("");
        });

/*
            if (protein[d.target.index].position !== protein[d.source.index].position) {
                return [numberWithCommas(d.source.value),
                        " trips from ",
                        protein[d.source.index].position,
                        " to ",
                        protein[d.target.index].position,
                        "\n",
                        numberWithCommas(d.target.value),
                        " trips from ",
                        protein[d.target.index].position,
                        " to ",
                        protein[d.source.index].position
                        ].join("");
                    //joining an array of many strings is faster than
                    //repeated calls to the '+' operator,
                    //and makes for neater code!
            }
            else { //source and target are the same
                return numberWithCommas(d.source.value)
                    + " trips started and ended in "
                    + protein[d.source.index].position;
            }
        });*/

    //handle exiting paths:
    chordPaths.exit().transition()
        .duration(250)
        .style("opacity", 1e-6)
        .remove();

    //update the path shape
    chordPaths
        .attr("opacity", hideZeros)
        .transition().duration(1000)
        .style("fill", function (d) {
            return chordFill(interactions[d.source.index][d.target.index]);
        })
        .attrTween("d", chordTween(last_layout))
        .transition().duration(100).attr("opacity", hideZeros)
    ;

    groupG.on("mouseover", function(d) {
        chordPaths.classed("fade", function (p) {
            return ((p.source.index != d.index) && (p.target.index != d.index));
        });
    });

    last_layout = layout; //save for next update

  });
}

function arcTween(oldLayout) {
    var oldGroups = {};
    if (oldLayout) {
        oldLayout.groups().forEach( function(groupData) {
            oldGroups[ groupData.index ] = groupData;
        });
    }

    return function (d, i) {
        var tween;
        var old = oldGroups[d.index];
        if (old) { //there's a matching old group
            tween = d3.interpolate(old, d);
        }
        else {
            //create a zero-width arc object
            var emptyArc = {startAngle: (d.startAngle + d.endAngle)*0.5,
                            endAngle: (d.startAngle + d.endAngle)*0.5};
            tween = d3.interpolate(emptyArc, d);
        }

        return function (t) {
            return arc( tween(t) );
        };
    };
}

function chordKey(data) {
    return (data.source.index < data.target.index) ?
        data.source.index  + "-" + data.target.index:
        data.target.index  + "-" + data.source.index;

    //create a key that will represent the relationship
    //between these two groups *regardless*
    //of which group is called 'source' and which 'target'
}
function chordTween(oldLayout) {
    //this function will be called once per update cycle

    //Create a key:value version of the old layout's chords array
    //so we can easily find the matching chord
    //(which may not have a matching index)

    var oldChords = {};

    if (oldLayout) {
        oldLayout.chords().forEach( function(chordData) {
            oldChords[ chordKey(chordData) ] = chordData;
        });
    }

    return function (d, i) {
        //this function will be called for each active chord

        var tween;
        var old = oldChords[ chordKey(d) ];
        if (old) {
            //old is not undefined, i.e.
            //there is a matching old chord value

            //check whether source and target have been switched:
            if (d.source.index != old.source.index ){
                //swap source and target to match the new data
                old = {
                    source: old.target,
                    target: old.source
                };
            }

            tween = d3.interpolate(old, d);
        }
        else {
            //create a zero-width chord object
            var emptyChord = {
                source: { startAngle: d.source.startAngle,
                         endAngle: d.source.startAngle},
                target: { startAngle: d.target.startAngle,
                         endAngle: d.target.startAngle}
            };
            tween = d3.interpolate( emptyChord, d );
        }

        return function (t) {
            //this function calculates the intermediary shapes
            return path(tween(t));
        };
    };
}

function disableAllButtons() {
    d3.selectAll("button").attr("disabled", "true");
}
function disableButton(buttonNode) {
    d3.selectAll("button").attr("disabled", function(d) {
        return this === buttonNode? "true": null;
    });
}
function enableTextSection(textNode) {
    var t = d3.select("pre"+textNode)[0][0].innerText;
    console.log(t);
    d3.select("#info").selectAll("p.descript")
        .transition().duration(200)
            .text(t)
            .style("margin-bottom", "50px")
            .style("font-size", "15px");
}

function callInitOnce() { if ( ! inited ) { inited = true; init(); } }
