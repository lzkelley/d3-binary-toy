/*
 *
 *
 *
 */
console.log("script.js");

/*    ========    SETTINGS    =========    */

var BINARY_MOTION = true;
var PARTICLE_MOTION = true;
var CONSUME = true;
var EVOLVE = false;
var EVOLUTION_INTERVAL = 20;
var THREE_D = false;

var NUM_PARTICLES = 400;

var svgSim = d3.select('#simContainer');
var svgPlot = d3.select('#plotContainer');

// These are strings
var simWidth = svgSim.style("width");
var simHeight = svgSim.style("height");
var plotWidth = svgPlot.style("width");
var plotHeight = svgPlot.style("height");

simWidth = parseFloat(simWidth.replace("px", ""));
simHeight = parseFloat(simHeight.replace("px", ""));
var simDepth = parseFloat(Math.max(simWidth, simHeight))

plotWidth = parseFloat(plotWidth.replace("px", ""));
plotHeight = parseFloat(plotHeight.replace("px", ""));

var VEL_SCALE = 0.02;
// var VEL_SCALE = 0.0;
var GRAVITY = 20.0
// var GRAVITY = 0.0;
var GRAVITY_SOFT_LENGTH = 20.0;
var PARTICLE_RADIUS = simWidth/200.0;
var NDIM = 2 + THREE_D;
var STAR_MASS = 0.01;
var STORE_TIME_INTERVAL_FACTOR = 1.0;
var XEDGE_INIT = 2.0;
var YEDGE_INIT = 2.0;

var LOG_SCALE = true;

if(NDIM == 2) {
    var bounds = [simWidth, simHeight];
} else {
    var bounds = [simWidth, simHeight, simDepth];
}

var maxSize = 50;
var sma = 100;
var period = 12*1000;
var sim_time = 0.0;
var store_time = 0.0;
var STAR_COLOR = "rgba(126, 126, 126, 0.5)";
var BH_COLOR = "#292929";
var BH_RADIUS = 40.0;
var padLR = 4;
var padTB = 4;
var BUTTON_PAD_X = 6;
var BUTTON_PAD_Y = 4;
var xedge_log_min = EVOLUTION_INTERVAL;
var yedge_log_min = 0.01;

var STATE_RUN = false;
var STATE_RESET = false;
var store_time_interval = EVOLUTION_INTERVAL * STORE_TIME_INTERVAL_FACTOR;

/*    ========    INITIALIZE OBJECTS    =========    */
// These are all set in `reset()`
var phaseInit, binary, particles, bhMassRatio, bhMassTotal, masses;

counterData = [
    {name: "Particles", num: NUM_PARTICLES, x: 0, y: 20},
    {name: "Consumed", num: 0, x: 0, y: 40},
    {name: "Created", num: 0, x: 0, y: 60}
]


// Storage for definitions, 'defs' tells SVG they are resources
var defs = svgSim.append("defs");

function bhSize(mass) {
    return mass*BH_RADIUS;
}

if(THREE_D){
    var pSizeScale = d3.scaleLinear()
    	.domain([0, simDepth])
    	.range([0.5*PARTICLE_RADIUS, 1.5*PARTICLE_RADIUS]);
    var pSize = function(pVec) { return pSizeScale(pVec[2]); }
} else {
    var pSize = function(pVec) { return PARTICLE_RADIUS; }
}

var gradientRadial = defs  // selectAll("radialGradient").data(binary).enter()
	.append("radialGradient")
	.attr("id", "radialGradient")
	.attr("cx", "30%")
	.attr("cy", "30%")
	.attr("r", "80%");

// Append the color stops
gradientRadial.append("stop")
	.attr("offset", "0%")
	.attr("stop-color", function(d) { return d3.rgb(BH_COLOR).brighter(2); });
gradientRadial.append("stop")
	.attr("offset", "50%")
	.attr("stop-color", function(d) { return BH_COLOR; });
gradientRadial.append("stop")
	.attr("offset",  "100%")
	.attr("stop-color", function(d) { return d3.rgb(BH_COLOR).darker(1.5); });

// == Integration Variables ==
var br = 0.0, bhRad = 0.0, soft = 0.0, distCubed = 0.0, pRad = 0.0;
var dr = new Array(NDIM);
var t0 = 0.0;
var t1 = 0.0;
var rad = 0.0, vals;
var lastTime = 0.0, edge = 0.0, newEdge = 0.0;

/*    ========    FUNCTIONS    =========    */

function format(num) {
    return num.toExponential(2);
}

function binaryPosition(bin){
    rad = sma * (bhMassTotal - bin.mass)/bhMassTotal
    vals = [
        0.5*simWidth + rad*Math.cos(2*Math.PI*bin.phase),
        0.5*simHeight + rad*Math.sin(2*Math.PI*bin.phase)
    ];
    if(THREE_D){
        vals.push(simDepth*0.5);
    }
    return vals;
}

function updateBinary(binary) {
    var selection = svgSim.selectAll(".bh")
    	.data(binary);

    // Remove old
    var olds = selection.exit();
    olds.remove();

    // for(var k = 0; k < 2; k++){
    //     console.log(k, binary[k].mass, binaryPosition(binary[k]));
    // }

    // Add new elements
    var news = selection.enter();
    news.append("circle")
        .attr("id", function(d){ return "bh-" + d.name; })
    	.attr("class", "bh")
    	.attr("r", function(d) { return bhSize(d.mass); })
    	.style("fill", function(d) { return "url(#radialGradient)"; })
        .attr("cx", function(d, i) { return binaryPosition(d)[0]; })
        .attr("cy", function(d, i) { return binaryPosition(d)[1]; });

    // Update elements
    // console.log("olds = ", selection.size());
    selection.attr("cx", function(d, i) { return binaryPosition(d)[0]; })
        .attr("cy", function(d, i) { return binaryPosition(d)[1]; })
        .attr("r", function(d) { return bhSize(d.mass); });
}

function updateParticles(particles) {
    // console.log("particles length = ", particles.length)
    var selection = svgSim.selectAll(".particleCircle")
        .data(particles);
    counterData[0].num = selection.size()
    // console.log("selection length = ", selection.size())

    // Update elements
    // console.log("selection: ", selection.size());
    selection.attr("cx", function(d, i) { return d[0]; })
        .attr("cy", function(d, i) { return d[1]; })
        .attr("r", function(d) { return pSize(d); })

    // Remove deleted elements
    var dels = selection.exit();
    // console.log("dels: ", dels.size());
    dels.remove();
    counterData[1].num += dels.size();

    // Add new elements
    var adds = selection.enter();
    // console.log("adds: ", adds.size());
    counterData[0].num += adds.size()
    adds.append("circle")
        .attr("id", function(d){ return "particleCircle-" + d.name; })
        .attr("class", "particleCircle")
        // .attr("r", PARTICLE_RADIUS)
        .attr("r", function(d) { return pSize(d); })
        .style("fill", STAR_COLOR)
        .style("opacity", 0.75)
        .attr("cx", function(d, i) { return d[0]; })
        .attr("cy", function(d, i) { return d[1]; })

};

function integrateBinary(binary, dt) {
    binary[0].phase += dt/period;
    binary[1].phase += dt/period;
    bhMassTotal = binary[0].mass + binary[1].mass;
}

function integrateParticles(particles, dt) {
    for(var ii = 0; ii < particles.length; ii++) {
        for(var jj = 0; jj < NDIM; jj++) {
            // if(ii == 0) console.log("pos =", particles[ii][jj]);
            // Update particle positions
            particles[ii][jj] += particles[ii][jj+NDIM] * dt * VEL_SCALE;
            // Update velocities
            particles[ii][jj+NDIM] += particles[ii][jj+2*NDIM] * dt * GRAVITY;
            // if(ii == 0) console.log("vel =", particles[ii][jj+NDIM]);
            // Check for out of bounds
            if(particles[ii][jj] < 0.0) { particles[ii][jj] += bounds[jj]; }
            if(particles[ii][jj] > bounds[jj]) { particles[ii][jj] -= bounds[jj]; }
        }

        // Update Accelerations
        for(var kk = 0; kk < 2; kk++) {
            br = binaryPosition(binary[kk]);
            bhRad = bhSize(binary[kk].mass)
            pRad = pSize(particles[ii])
            // If `CONSUME` then we dont need softening
            soft = bhRad * !CONSUME;

            distCubed = Math.pow(soft, 2);

            for(var jj = 0; jj < NDIM; jj++) {
                dr[jj] = particles[ii][jj] - br[jj];
                distCubed += Math.pow(dr[jj], 2.0);
            }
            distCubed = Math.pow(distCubed, 1.5);

            if(CONSUME && (Math.pow(bhRad + pRad, 3.0) >= distCubed)){
                particles.splice(ii, 1);
                binary[kk].mass += STAR_MASS;
                ii--;
                break;
            } else {
                for(var jj = 0; jj < NDIM; jj++) {
                    particles[ii][jj + 2*NDIM] = - dr[jj] * binary[kk].mass/distCubed;
                }
            }
        }

        // console.log(ii, particles[ii]);
        // if(isNaN(particles[ii][0])) {
        //     console.log("STOP");
        //     interval.stop();
        //     break;
        // }

    }
}

function store(tt) {
    masses.push([tt/1000.0, binary[0].mass, binary[1].mass]);
}

function evolve(tt) {
    if(!STATE_RUN) {
        return;
    }

    t0 = t1;
    t1 = tt;
    dt = (t1 - t0) * (t0 > 0 && t1 > 0);
    sim_time += dt;
    store_time += dt;

    if(store_time > store_time_interval){
        store_time = 0.0;
        store(sim_time);
        // console.log("Length = ", masses.length);
    }

    if(PARTICLE_MOTION) {
        integrateParticles(particles, dt);
        updateParticles(particles);
    }

    if(BINARY_MOTION) {
        integrateBinary(binary, dt);
        updateBinary(binary);
    }

    updatePlots();
    updateCounters();
}

function reset() {
    console.log("RESET!");

    bhMassTotal = 1.0;
    bhMassRatio = Math.random();
    m1 = bhMassTotal/(1.0 + bhMassRatio);
    m2 = bhMassTotal - m1;
    console.log("m1 = ", format(m1), "; m2 = ",
                format(m2), "; mu = ", format(m2/m1), "; M = ", format(m1+m2));

    phaseInit = Math.random();
    binary = [
    	{name: "a", phase: phaseInit, mass: m1},
    	{name: "b", phase: phaseInit + 0.5, mass: m2}
    ];

    particles = d3.range(NUM_PARTICLES).map(function(i) {
        var vals = new Array(3*NDIM);
        for(var ii = 0; ii < NDIM; ii++) {
            // Positions
            vals[ii] = bounds[ii] * Math.random();
            // Velocities
            vals[ii + NDIM] = 2*(Math.random() - 0.5);
            // Accelerations
            vals[ii + NDIM*2] = 0.0;
        }
        return vals;
    });

    counterData[0].num = NUM_PARTICLES;
    counterData[1].num = 0;
    counterData[2].num = 0;

    sim_time = 0.0;
    store_time = 0.0;
    masses = [[0.0, m1, m2]];

    updateCounters();
    updateParticles(particles);
    updateBinary(binary);
    updatePlots(true);
}

function initCounters() {
    var texts = svgSim.selectAll("text.counters")
        .data(counterData).enter();

    texts.append("rect")
        .attr("id", function(d){ return d.name; })
        .attr("class", "counters")
        .attr("x", function (d, i) { return d.x; })
        .attr("y", function (d, i) { return d.y - 20; });

    texts.append("text")
        .attr("id", function(d){ return d.name; })
        .attr("class", "counters")
        .attr("text-anchor", "left")
        .attr("x", function (d, i) { return d.x; })
        .attr("y", function (d, i) { return d.y; });
}

function updateCounters() {

    texts = svgSim.selectAll("text.counters");

    texts.text(function(d) { return d.name + ": " + d.num; });

    // get bounding box of text field and store it in texts array
    texts.each(function(d, i) { d.bb = this.getBBox(); });

    svgSim.selectAll("rect.counters")
        .attr("x", function(d) { return d.x - padLR/2; })
        .attr("y", function(d) { return d.y + padTB/2 - 20;  })
        .attr("simWidth", function(d) { return d.bb.width + padLR; })
        .attr("simHeight", function(d) { return d.bb.height + padTB; });
}

function initButtons() {
    // On/Off Button
    svgSim.append("rect")
        .attr("class", "button")
        .attr("id", "stateRect")
        .on("click", toggleState);

    svgSim.append("text")
        .attr("class", "button")
        .attr("id", "stateText")
        .attr("text-anchor", "end")
        .attr("x", "98%")
        .attr("y", "2%")
        .text("Start")
        .on("click", toggleState);

    // Reset Button
    svgSim.append("rect")
        .attr("class", "button")
        .attr("id", "resetRect")
        .on("click", toggleState);

    svgSim.append("text")
        .attr("class", "button")
        .attr("id", "resetText")
        .attr("text-anchor", "end")
        .attr("x", "98%")
        .attr("y", "4%")
        .text("Reset")
        .on("click", reset);

}

function toggleState(d, i) {
    STATE_RUN = !STATE_RUN
    updateButtons();
    t0 = t1 = 0.0;
}

function updateButtons() {

    // == State Button == //
    stateText = svgSim.selectAll("#stateText");
    textBB = stateText.node().getBBox();
    var yy = 0.02*simHeight + textBB.height;
    stateText.attr("y", yy);

    stateText.text(function (d) {
        if (!STATE_RUN) {
            return "Start";
        } else {
            return "Stop";
        }
    });

    // Have to update the BBox
    textBB = stateText.node().getBBox();
    stateRect = svgSim.selectAll("#stateRect")
        .style("x", textBB.x - BUTTON_PAD_X)
        .style("y", textBB.y - BUTTON_PAD_Y)
        .style("width", textBB.width + 2*BUTTON_PAD_X)
        .style("height", textBB.height + 2*BUTTON_PAD_Y);

    // == Reset Button == //
    rectBB = stateRect.node().getBBox();
    resetText = svgSim.selectAll("#resetText");
    // textBB = stateText.node().getBBox();
    var yy = 0.02*simHeight + textBB.height + rectBB.height;
    resetText.attr("y", yy);

    // Have to update the BBox
    textBB = resetText.node().getBBox();
    resetRect = svgSim.selectAll("#resetRect")
        .style("x", textBB.x - BUTTON_PAD_X)
        .style("y", textBB.y - BUTTON_PAD_Y)
        .style("width", textBB.width + 2*BUTTON_PAD_X)
        .style("height", textBB.height + 2*BUTTON_PAD_Y);

}

function initPlots() {
    var pathStringM1 = lineGenM1(masses);
    var pathStringM2 = lineGenM2(masses);

    svgPlot.append("path")
        .attr("class", "line")
        .attr("id", "m1")
        .attr("transform", "translate(50, 10)")
    	.attr('d', pathStringM1);

    svgPlot.append("path")
        .attr("class", "line")
        .attr("id", "m2")
        .attr("transform", "translate(50, 10)")
    	.attr('d', pathStringM2);

}

function updatePlots(reset=false) {

    // == Update Axes == //

    // Update xaxes
    lastTime = masses.slice(-1)[0][0];
    edge = xscale.domain()[1];
    if(lastTime > 0.8*edge || reset){
        if (reset) {
            newEdge = XEDGE_INIT;
        } else {
            newEdge = 2.0 * edge;
        }
        xscale.domain([0.0, newEdge]);

        svgPlot.select("#xaxis")
            .transition().duration(500).ease(d3.easeSinInOut)
            .call(xAxis);
    }

    // Update yaxes
    lastMass = masses.slice(-1)[0];
    console.log(lastMass);
    lastMass = Math.max(lastMass[1], lastMass[2]);
    console.log(lastMass);
    edge = yscale.domain()[1];
    if(lastMass > 0.8*edge || reset){
        if (reset) {
            newEdge = YEDGE_INIT;
        } else {
            newEdge = 1.5 * edge;
        }
        yscale.domain([0.0, newEdge]);

        svgPlot.select("#yaxis")
            .transition().duration(500).ease(d3.easeSinInOut)
            .call(yAxis);
    }

    // == Update Lines == //

    // SVG string specification for the line
    var pathStringM1 = lineGenM1(masses);
    var pathStringM2 = lineGenM2(masses);

    svgPlot.select("#m1")
    	.attr('d', pathStringM1);

    svgPlot.select("#m2")
    	.attr('d', pathStringM2);
}

/*   =========    PLOTS        =========   */

var xscale = d3.scaleLinear()
               .domain([0.0, XEDGE_INIT])
               .range([0.0, plotWidth - 100]);

var yscale = d3.scaleLinear()
               .domain([0, YEDGE_INIT])
               .range([plotHeight/2, 0]);

var xAxis = d3.axisBottom().scale(xscale);

var yAxis = d3.axisLeft().ticks(5).scale(yscale);

svgPlot.append("g")
    .attr("class", "axis")
    .attr("id", "yaxis")
    .attr("transform", "translate(50, 10)")
    .call(yAxis);

var xAxisTranslate = plotHeight/2 + 10;

svgPlot.append("g")
    .attr("class", "axis")
    .attr("id", "xaxis")
    .attr("transform", "translate(50, " + xAxisTranslate  +")")
    .call(xAxis);

var lineGenM1 = d3.line()
    .x(function(d) { return xscale(d[0]); })
    .y(function(d) { return yscale(d[1]); });

var lineGenM2 = d3.line()
    .x(function(d) { return xscale(d[0]); })
    .y(function(d) { return yscale(d[2]); });


/*    ========    SIMULATION    =========    */

reset();

// console.log("Masses = ", masses);

initPlots();

// Call this after particles and binaries so that it's on top.
initCounters();
updateCounters();

initButtons();
updateButtons();

var interval = d3.interval(evolve, EVOLUTION_INTERVAL);
