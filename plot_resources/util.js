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

Util.transpose = function(mat) {
  let M = mat.length;
  let N = mat[0].length;
  let T = [];

  for(let n = 0; n < N; n++) {
    T.push([]);
    for(let m = 0; m < M; m++) {
      T[n].push(mat[m][n]);
    }
  }
  return T;
}

Util.find_clonal = function(parents) {
  let clonal = [];
  for(var idx = 0; idx < parents.length; idx++) {
    if(parents[idx] === 0) {
      clonal.push(idx + 1);
    }
  }
  return new Set(clonal);
}

Util.calc_ccf = function(phi, parents) {
  let S = phi[0].length;
  let clonal = Util.find_clonal(parents);
  let clonal_phi = [];

  for(let sidx = 0; sidx < S; sidx++) {
    clonal_phi.push(0);
    clonal.forEach(cidx => clonal_phi[sidx] += phi[cidx][sidx]);
    console.assert(clonal_phi[sidx] <= 1);
  }

  let ccf = phi.slice(1).map(P_k => {
    return P_k.map((P_ks, sidx) => P_ks / clonal_phi[sidx]);
  });

  return ccf;
}


function ColourAssigner() {
}

ColourAssigner._assign_from_palette = function(num_colours, scale_name) {
  let scales = {
    tableau10: d3.schemeTableau10,
    dark24: ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038'],
    alphabet: ['#2E91E5', '#E15F99', '#1CA71C', '#FB0D0D', '#DA16FF', '#222A2A', '#B68100', '#750D86', '#EB663B', '#511CFB', '#00A08B', '#FB00D1', '#FC0080', '#B2828D', '#6C7C32', '#778AAE', '#862A16', '#A777F1', '#620042', '#1616A7', '#DA60CA', '#6C4516', '#0D2A63', '#AF0038'],
  };
  let scale = scale_name ? scales[scale_name] : scales.tableau10;

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

ColourAssigner._assign_from_husl = function(N) {
  let colours = [];
  let first_hue = 0.01;
  for(let idx = 0; idx < N; idx++) {
    let hue = 360*(first_hue + ((idx+1)/N)*(1 - first_hue));
    // Docs: https://github.com/hsluv/hsluv/tree/master/javascript
    colours.push(hsluv.hsluvToHex([hue, 66, 56]));
  }
  shuffle(colours);
  return colours;
}

ColourAssigner.assign_colours = function(N, use_husl=false) {
  if(use_husl) {
    return ColourAssigner._assign_from_husl(N);
  } else {
    return ColourAssigner._assign_from_palette(N);
  }
}

function shuffle(array) {
  for (let i = array.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [array[i], array[j]] = [array[j], array[i]];
  }
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
