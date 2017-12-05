function TreePlotter() { }

// Following two functions from https://stackoverflow.com/a/18473154.
function polarToCartesian(centerX, centerY, radius, angleInDegrees) {
  var angleInRadians = (angleInDegrees-90) * Math.PI / 180.0;

  return {
    x: centerX + (radius * Math.cos(angleInRadians)),
    y: centerY + (radius * Math.sin(angleInRadians))
  };
}

function describeArc(x, y, radius, startAngle, endAngle){
  var start = polarToCartesian(x, y, radius, endAngle);
  var end = polarToCartesian(x, y, radius, startAngle);

  var largeArcFlag = endAngle - startAngle <= 180 ? "0" : "1";

  var d = [
    "M", start.x, start.y,
    "A", radius, radius, 0, largeArcFlag, 0, end.x, end.y
  ].join(" ");

  return d;
}

TreePlotter.prototype._calculate_max_depth = function(root) {
  var _calc_max_depth = function(node) {
    if(!node.hasOwnProperty('children')) {
      return 0;
    }
    var max_depth = 0;
    node.children.forEach(function(child) {
      var nd = _calc_max_depth(child);
      if(nd > max_depth)
        max_depth = nd;
    });
    return 1 + max_depth;
  };
  return _calc_max_depth(root);
}

TreePlotter.prototype._label_node = function(node_id) {
  // If we're examining a patient+xeno tree, label the node with the node_id.
  // Otherwise, we're presumably examining a patient-only tree, so label it
  // with a letter.
  if(window.location.toString().indexOf('pairwise.xeno') > -1) {
    return node_id;
  }
  // Restrict to valid alphabet range.
  if(node_id < 1 || node_id > 26) {
    return node_id;
  }
  var letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  return letters[node_id - 1];
}

TreePlotter.prototype._draw_tree = function(root, container, num_pops) {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  var max_depth = this._calculate_max_depth(root);
  var m = [10, horiz_padding, 10, horiz_padding],
      w = 120*max_depth - m[1] - m[3],
      h = 600 - m[0] - m[2],
      i = 0;
  var colours = ColourAssigner.assign_colours(num_pops);

  var diag_colour    = '#428bca';
  var relapse_colour = '#ca4242';
  var node_bgcolour  = '#ffffff';

  // Compute the new tree layout.
  var tree = d3.tree().size([h, w]);
  root = tree(d3.hierarchy(root));
  root.descendants().sort(function(a, b) {
    return d3.ascending(parseInt(a.data.name, 10), parseInt(b.data.name, 10));
  });
  var svg = container.append('svg:svg')
      .attr('width', w + m[1] + m[3])
      .attr('height', h + m[0] + m[2]);
  var vis = svg.append('svg:g')
      .attr('transform', 'translate(' + m[3] + ',' + m[0] + ')');

  // Update the nodes.
  var node = vis.selectAll('g.node')
      .data(root.descendants(), function(d) { return d.data.name; });

  var self = this;
  // Enter any new nodes at the parent's previous position.
  var nodeEnter = node.enter().append('svg:g');
  nodeEnter.attr('class', 'population')
    .attr('fill', function(d, i) { return colours[d.data.name]; })
    .attr('transform', function(d) { return 'translate(' + d.y + ',' + d.x + ')'; });

  nodeEnter.append('svg:circle')
    .attr('class', 'outline')
    .attr('r', function(d) { return d.data.radius; });
  nodeEnter.append('svg:path')
    .attr('class', 'left_half')
    .attr('d', function(d) { return describeArc(0, 0, d.data.radius, 180, 360); })
    .attr('fill', function(d) { return Util.rgba2hex(diag_colour, d.data.opacities.left, node_bgcolour); });
  nodeEnter.append('svg:path')
    .attr('class', 'right_half')
    .attr('d', function(d) { return describeArc(0, 0, d.data.radius, 0, 180); })
    .attr('fill', function(d) { return Util.rgba2hex(relapse_colour, d.data.opacities.right, node_bgcolour); });
  nodeEnter.append('svg:path')
    .attr('class', 'divider')
    .attr('stroke', '#aaa')
    .attr('d', function(d) { return 'M 0 -' + d.data.radius + ' V ' + d.data.radius; });

  nodeEnter.append('svg:text')
      .attr('font-size', '30')
      .attr('dominant-baseline', 'central')
      .attr('text-anchor', 'middle')
      .text(function(d) { return self._label_node(d.data.name); });

  // Update the links.
  var link = vis.selectAll('path.link')
      .data(root.links(), function(d) { return d.target.data.name; })
      .attr('stroke-width', '1.5px')

  // Enter any new links at the parent's previous position.
  link.enter().insert('svg:path', 'g')
    .attr('class', 'link')
    .attr('stroke', '#aaa')
    .attr('d', d3.linkHorizontal().x(function(d) {
      return d.y;
    }).y(function(d) {
      return d.x;
    }));
}

TreePlotter.prototype._find_max_ssms = function(populations) {
  var max_ssms = 0;
  for(var pop_id in populations) {
    var pop = populations[pop_id];
    if(pop.num_ssms > max_ssms)
      max_ssms = pop.num_ssms;
  }
  return max_ssms;
}

TreePlotter.prototype._generate_tree_struct = function(sampnames, adjlist, pops, root_id, left_sample, right_sample) {
  var max_ssms = this._find_max_ssms(pops);

  var left_index = sampnames.indexOf(left_sample);
  var right_index = sampnames.indexOf(right_sample);
  if(left_index === -1 || right_index === -1) {
    throw "Can't find " + left_sample + " or " + right_sample + " samples in " + sampnames;
  }

  var _add_node = function(node_id, struct) {
    struct.name = node_id;

    var num_ssms = pops[node_id]['num_ssms'];
    struct.radius = TreeUtil.calc_radius(num_ssms /  max_ssms);
    struct.opacities = {
      'left':   pops[node_id]['cellular_prevalence'][left_index],
      'right':  pops[node_id]['cellular_prevalence'][right_index],
    };

    if(typeof adjlist[node_id] === 'undefined') {
      return;
    }
    struct.children = [];
    adjlist[node_id].forEach(function(child_id) {
      var child = {};
      struct.children.push(child);
      _add_node(child_id, child);
    });
  };

  var root = {};
  _add_node(root_id, root);
  return root;
}

TreePlotter.prototype.plot = function(summ_path, tidx, tname, left_sample, right_sample, container) {
  container = d3.select(container).append('div');
  container.append('h2').text(tname + ' tidx=' + tidx + ' left=' + left_sample + ' right=' + right_sample);

  var self = this;
  d3.json(summ_path, function(summary) {
    var pops = summary.trees[tidx].populations;
    var struct = summary.trees[tidx].structure;
    var root = summary.trees[tidx].root;
    if(struct[root].length !== 1) {
      throw "Unexpected children from root: " + struct[root];
    }
    var clonal = struct[root][0];

    var root = self._generate_tree_struct(summary.params.samples, struct, pops, clonal, left_sample, right_sample);
    self._draw_tree(root, container, Object.keys(pops).length);
  });
}

function TreeUtil() {
}

TreeUtil.calc_radius = function(scale) {
  var min_area = 700, max_area = 8000;
  // scale must be in [0, 1].
  var area = min_area + scale*(max_area - min_area);
  return Math.sqrt(area / Math.PI);
}

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

Util.rgba2hex = function(base_colour, alpha, bgcolour) {
  base_colour = Util.hex2rgb(base_colour);
  bgcolour = Util.hex2rgb(bgcolour);

  var rgb = {};
  Object.keys(base_colour).forEach(function(channel) {
    rgb[channel] = Math.round(alpha*base_colour[channel] + (1 - alpha)*bgcolour[channel]);
  });
  var hex = '#' + rgb.r.toString(16) + rgb.g.toString(16) + rgb.b.toString(16);
  console.log([base_colour, alpha, bgcolour, rgb, hex]);
  return hex;
}

function PhiMatrix() {
}

PhiMatrix.prototype._calc_ccf = function(phi) {
  var ccf = [];
  var clonalidx = 1;

  for(var rowidx = clonalidx; rowidx < phi.length; rowidx++) {
    ccf[rowidx - 1] = [];
    for(var sampidx = 0; sampidx < phi[rowidx].length; sampidx++) {
      ccf[rowidx - 1].push(phi[rowidx][sampidx] / phi[clonalidx][sampidx]);
    }
  }

  return ccf;
}

PhiMatrix.prototype.plot = function(phi_path, container) {
  var self = this;
  d3.json(phi_path, function(phi_data) {
    var ccf = self._calc_ccf(phi_data['phi']);
    var sampnames = phi_data['samples'];

    var num_rows = ccf.length;
    var num_cols = ccf[0].length;
    var cell_size = 50;
    var row_label_width = 100;
    var col_label_height = 100;
    var label_padding = 10;

    var colours = ColourAssigner.assign_colours(phi_data['phi'].length);

    var svg = d3.select(container).html('').append('svg:svg')
      .attr('width', row_label_width + num_cols * cell_size)
      .attr('height', col_label_height + num_rows * cell_size);
    svg.append('svg:g')
      .attr('transform', function(d, i) { return 'translate(' + (row_label_width + 0.5 * cell_size) + ',' + (col_label_height - label_padding) + ')'; })
      .selectAll('text')
      .data(sampnames)
      .enter()
      .append('svg:text')
      .attr('transform', function(d, i) { return 'translate(' + i * cell_size + ',0) rotate(270)'; })
      .attr('x', 0)
      .attr('y', 0)
      .text(function(d, i) { return d; });
    var rows = svg.selectAll('g.phis')
      .data(ccf)
      .enter()
      .append('svg:g')
      .attr('class', 'phis')
      .attr('fill', function(d, i) { return colours[i + 1]; })
      .attr('transform', function(d, i) { return 'translate(' + row_label_width + ',' + (col_label_height + (i * cell_size)) + ')'; });
    rows.append('svg:text')
      .attr('x', -label_padding)
      .attr('y', 0.5 * cell_size)
      .attr('dominant-baseline', 'middle')
      .attr('text-anchor', 'end')
      .text(function(d, i) { return 'Population ' + (i + 1); });
    rows.selectAll('rect')
      .data(function(d) { return d; })
      .enter()
      .append('svg:rect')
      .attr('width', cell_size)
      .attr('height', function(d, i) { return d*cell_size; })
      .attr('x', function(d, i) { return i * cell_size; })
      .attr('y', function(d, i) { return (1 - d)*cell_size; })
      .attr('fill-opacity', function(d) { return 1.0; });
  });
}

function ColourAssigner() {
}

ColourAssigner.assign_colours = function(num_colours) {
  var scale = d3.schemeDark2;

  var colours = [];
  for(var cidx = 0; colours.length < num_colours; cidx++) {
    if(cidx === scale.length) {
      cidx = 0;
    }
    colours.push(scale[cidx]);
  }
  return colours;
}