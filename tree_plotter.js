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

TreePlotter.prototype._draw_tree = function(root, container) {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  var m = [10, horiz_padding, 10, horiz_padding],
      w = 800 - m[1] - m[3],
      h = 600 - m[0] - m[2],
      i = 0;

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

  // Enter any new nodes at the parent's previous position.
  var nodeEnter = node.enter().append('svg:g');
  nodeEnter.attr('class', 'population')
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

TreePlotter.prototype._draw = function(populations, structure, root_id, params) {
  if(params !== undefined) { }

  var root = this._generate_tree_struct(structure, populations, root_id);
  this._draw_tree(root);
}

TreePlotter.prototype.plot = function(summ_path, tidx, container) {
  var self = this;
  d3.json(summ_path, function(summary) {
    var root = self._generate_tree_struct(summary.params.samples, summary.trees[tidx].structure, summary.trees[tidx].populations, summary.trees[tidx].root);
    self._draw_tree(root, container);
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
