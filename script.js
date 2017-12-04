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
var EVOLUTION_INTERVAL = 10;
var THREE_D = false;

var NUM_PARTICLES = 800;

var svg = d3.select('#simContainer');

// These are strings
var width = svg.style("width");
var height = svg.style("height");

width = parseFloat(width.replace("px", ""));
height = parseFloat(height.replace("px", ""));
var depth = parseFloat(Math.max(width, height))

var VEL_SCALE = 0.02;
// var VEL_SCALE = 0.0;
var GRAVITY = 20.0
// var GRAVITY = 0.0;
var GRAVITY_SOFT_LENGTH = 20.0;
var PARTICLE_RADIUS = width/200.0;
var NDIM = 2 + THREE_D;
var STAR_MASS = 0.001;

if(NDIM == 2) {
    var bounds = [width, height];
} else {
    var bounds = [width, height, depth];
}

var maxSize = 50;
var sma = 100;
var period = 12*1000;
var time = 0.0;
var BH_COLOR = "#292929";
var BH_RADIUS = 40.0;
var padLR = 4;
var padTB = 4;
var BUTTON_PAD_X = 6;
var BUTTON_PAD_Y = 4;

var STATE_RUN = false;
var STATE_RESET = false;

/*    ========    INITIALIZE OBJECTS    =========    */
// These are all set in `reset()`
var phaseInit, binary, particles, bhMassRatio, bhMassTotal;

counterData = [
    {name: "Particles", num: NUM_PARTICLES, x: 0, y: 20},
    {name: "Consumed", num: 0, x: 0, y: 40},
    {name: "Created", num: 0, x: 0, y: 60}
]


// Storage for definitions, 'defs' tells SVG they are resources
var defs = svg.append("defs");

function bhSize(mass) {
    return mass*BH_RADIUS;
}

if(THREE_D){
    var pSizeScale = d3.scaleLinear()
    	.domain([0, depth])
    	.range([0.5*PARTICLE_RADIUS, 1.5*PARTICLE_RADIUS]);
    var pSize = function(pVec) { return pSizeScale(pVec[2]); }
} else {
    var pSize = function(pVec) { return PARTICLE_RADIUS; }
}

var gradientRadial = defs  //selectAll("radialGradient").data(binary).enter()
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

/*    ========    FUNCTIONS    =========    */

function format(num) {
    return num.toExponential(2);
}

function binaryPosition(bin){
    rad = sma * (bhMassTotal - bin.mass)/bhMassTotal
    vals = [
        0.5*width + rad*Math.cos(2*Math.PI*bin.phase),
        0.5*height + rad*Math.sin(2*Math.PI*bin.phase)
    ];
    if(THREE_D){
        vals.push(depth*0.5);
    }
    return vals;
}

function updateBinary(binary) {
    var selection = svg.selectAll(".bh")
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
    var selection = svg.selectAll(".particleCircle")
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
        .style("fill", "red")
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

function evolve(tt) {
    if(!STATE_RUN) {
        return;
    }

    t0 = t1;
    t1 = tt;
    dt = (t1 - t0) * (t0 > 0 && t1 > 0);
    time += dt;

    if(PARTICLE_MOTION) {
        integrateParticles(particles, dt);
        updateParticles(particles);
    }

    if(BINARY_MOTION) {
        integrateBinary(binary, dt);
        updateBinary(binary);
    }

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

    updateCounters();
    updateParticles(particles);
    updateBinary(binary);

}

function initCounters() {
    var texts = svg.selectAll("text.counters")
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

    texts = svg.selectAll("text.counters");

    texts.text(function(d) { return d.name + ": " + d.num; });

    // get bounding box of text field and store it in texts array
    texts.each(function(d, i) { d.bb = this.getBBox(); });

    svg.selectAll("rect.counters")
        .attr("x", function(d) { return d.x - padLR/2; })
        .attr("y", function(d) { return d.y + padTB/2 - 20;  })
        .attr("width", function(d) { return d.bb.width + padLR; })
        .attr("height", function(d) { return d.bb.height + padTB; });
}

function initButtons() {
    // On/Off Button
    svg.append("rect")
        .attr("class", "button")
        .attr("id", "stateRect")
        .on("click", toggleState);

    svg.append("text")
        .attr("class", "button")
        .attr("id", "stateText")
        .attr("text-anchor", "end")
        .attr("x", "98%")
        .attr("y", "2%")
        .text("Start")
        .on("click", toggleState);

    // Reset Button
    svg.append("rect")
        .attr("class", "button")
        .attr("id", "resetRect")
        .on("click", toggleState);

    svg.append("text")
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
    stateText = svg.selectAll("#stateText");
    textBB = stateText.node().getBBox();
    stateText.attr("y", 0.02*height + textBB.height);

    stateText.text(function (d) {
        if (!STATE_RUN) {
            return "Start";
        } else {
            return "Stop";
        }
    });

    // Have to update the BBox
    textBB = stateText.node().getBBox();
    stateRect = svg.selectAll("#stateRect")
        .style("x", textBB.x - BUTTON_PAD_X)
        .style("y", textBB.y - BUTTON_PAD_Y)
        .style("width", textBB.width + 2*BUTTON_PAD_X)
        .style("height", textBB.height + 2*BUTTON_PAD_Y);

    // == Reset Button == //
    rectBB = stateRect.node().getBBox();
    resetText = svg.selectAll("#resetText");
    textBB = stateText.node().getBBox();
    resetText.attr("y", 0.02*height + textBB.height + rectBB.height);

    // Have to update the BBox
    textBB = resetText.node().getBBox();
    resetRect = svg.selectAll("#resetRect")
        .style("x", textBB.x - BUTTON_PAD_X)
        .style("y", textBB.y - BUTTON_PAD_Y)
        .style("width", textBB.width + 2*BUTTON_PAD_X)
        .style("height", textBB.height + 2*BUTTON_PAD_Y);

}

/*    ========    SIMULATION    =========    */

reset();

// Call this after particles and binaries so that it's on top.
initCounters();
updateCounters();

initButtons();
updateButtons();

var interval = d3.interval(evolve, EVOLUTION_INTERVAL);
