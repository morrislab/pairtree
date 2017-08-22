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

  // Compute the new tree layout.
  var tree = d3.tree().size([h, w]);
  root = tree(d3.hierarchy(root));
  root.descendants().sort(function(a, b) {
    return d3.ascending(parseInt(a.data.name, 10), parseInt(b.data.name, 10));
  });
  var svg = d3.select(container).html('').append('svg:svg')
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
    .attr('class', 'diagnosis')
    .attr('d', function(d) { return describeArc(0, 0, d.data.radius, 180, 360); })
    .attr('opacity', function (d) { return d.data.patient_cell_prevs.diagnosis; });
  nodeEnter.append('svg:path')
    .attr('class', 'relapse')
    .attr('d', function(d) { return describeArc(0, 0, d.data.radius, 0, 180); })
    .attr('opacity', function (d) { return d.data.patient_cell_prevs.relapse; });
  nodeEnter.append('svg:path')
    .attr('class', 'divider')
    .attr('stroke', '#aaa')
    .attr('d', function(d) { return 'M 0 -' + d.data.radius + ' V ' + d.data.radius; });

  nodeEnter.append('svg:text')
      .attr('font-size', '30')
      .attr('dominant-baseline', 'central')
      .attr('text-anchor', 'middle')
      .text(function(d) { return d.data.name; });

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

TreePlotter.prototype._generate_tree_struct = function(sampnames, adjlist, pops, root_id) {
  var max_ssms = this._find_max_ssms(pops);

  var diag_index = sampnames.indexOf('D');
  var relapse_index = sampnames.indexOf('R1');
  if(diag_index === -1 || relapse_index === -1) {
    throw "Can't find diagnosis or relapse sample in " + sampnames;
  }

  var _add_node = function(node_id, struct) {
    struct.name = node_id;

    var num_ssms = pops[node_id]['num_ssms'];
    struct.radius = TreeUtil.calc_radius(num_ssms /  max_ssms);
    struct.patient_cell_prevs = {
      'diagnosis': pops[node_id]['cellular_prevalence'][diag_index],
      'relapse':   pops[node_id]['cellular_prevalence'][relapse_index],
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

TreePlotter.prototype.plot = function(summ_path, tidx, container) {
  var self = this;
  d3.json(summ_path, function(summary) {
    var root = self._generate_tree_struct(summary.params.samples, summary.trees[tidx].structure, summary.trees[tidx].populations, summary.trees[tidx].root);
    self._draw_tree(root, container, Object.keys(summary.trees[tidx].populations).length);
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
      .attr('height', function(d, i) { return cell_size; return d*cell_size; })
      .attr('x', function(d, i) { return i * cell_size; })
      .attr('y', function(d, i) { return 0; return 0.5*(1 - d)*cell_size; })
      .attr('fill-opacity', function(d) { return d; });
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
