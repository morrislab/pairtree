function Util() {
}

Util.sort_ints = function(arr) {
  var converted = arr.map(function(a) {
    return parseInt(a, 10);
  });
  return converted.sort();
}

Util.hex2rgb = function(colour) {
  return {
    r: parseInt(colour.substring(1, 3), 16),
    g: parseInt(colour.substring(3, 5), 16),
    b: parseInt(colour.substring(5, 7), 16)
  };
}

Util.blend_colours = function(base_colour, alpha, bgcolour) {
  base_colour = Util.hex2rgb(base_colour);
  bgcolour = Util.hex2rgb(bgcolour);

  var rgb = {};
  Object.keys(base_colour).forEach(function(channel) {
    rgb[channel] = Math.round(alpha*base_colour[channel] + (1 - alpha)*bgcolour[channel]);
    rgb[channel] = rgb[channel].toString(16).padStart(2, "0");
  });
  var hex = '#' + rgb.r + rgb.g + rgb.b;
  return hex;
}

function ColourAssigner() {
}

ColourAssigner.assign_colours = function(num_colours) {
  var scale = d3.schemeTableau10;

  var colours = [];
  // Start at cidx=1 rather than cidx=0 to keep colouring of phi matrix
  // consistent with what it was in past versions.
  for(var cidx = 1; colours.length < num_colours; cidx++) {
    if(cidx === scale.length) {
      cidx = 0;
    }
    colours.push(scale[cidx]);
  }
  return colours;
}

function resize_svg(elems) {
  elems.nodes().forEach(function(svg) {
    var box = svg.getBBox();
    // Padding is necessary, as `getBBox` doesn't account for borders around
    // nodes -- these would otherwise be clipped off.
    var padding = 6;
    var viewbox = [box.x - 0.5*padding, box.y - 0.5*padding, box.width + padding, box.height + padding];
    svg.setAttribute('viewBox', viewbox.join(' '));
    svg.setAttribute('width', viewbox[2]);
    svg.setAttribute('height', viewbox[3]);
  });
}
