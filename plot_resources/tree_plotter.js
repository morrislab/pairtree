// 'letters' or 'digits'
var LABEL_TYPE = 'digits';
// 'circle' or 'square'
var NODE_TYPE = 'circle';

// Following two functions from https://stackoverflow.com/a/18473154.
function polarToCartesian(centerX, centerY, radius, angleInDegrees) {
  // Add 90 degrees to rotate coordinate system so that arcs start at the
  // bottom of the circle and rotate clockwise.
  var angleInRadians = (angleInDegrees + 90) * Math.PI/180.0;
  return {
    x: centerX + (radius * Math.cos(angleInRadians)),
    y: centerY + (radius * Math.sin(angleInRadians))
  };
}

function describeArc(x, y, radius, startAngle, endAngle) {
  var start = polarToCartesian(x, y, radius, startAngle);
  var end = polarToCartesian(x, y, radius, endAngle);
  var largeArcFlag = endAngle - startAngle <= 180 ? "0" : "1";
  var d = [
    "M", start.x, start.y,
    "A", radius, radius, 0, largeArcFlag, 1, end.x, end.y,
    "L", x, y,
    "L", start.x, start.y,
  ].join(" ");

  return d;
}

function describeRect(x, y, width, height) {
  var d = [
    "M", x, y,
    "H", x + width,
    "V", y + height,
    "H", x,
    "Z"
  ].join(" ");
  return d;
}

function TreePlotter() {
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
  if(LABEL_TYPE === 'digits') {
    return node_id;
  }
  // Restrict to valid alphabet range.
  if(node_id < 1 || node_id > 26) {
    return node_id;
  }
  var letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
  return letters[node_id - 1];
}

TreePlotter.prototype._draw_tree = function(root, container, sampnames, samp_colours, pop_colours, bg_colour) {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  var max_depth = this._calculate_max_depth(root);
  var m = [10, horiz_padding, 10, horiz_padding],
      w = Math.max(150, 120*max_depth - m[1] - m[3]),
      h = 430 - m[0] - m[2],
      i = 0;

  // Compute the new tree layout.
  var tree = d3.tree().size([h, w]);
  root = tree(d3.hierarchy(root));
  root.descendants().sort(function(a, b) {
    return d3.ascending(parseInt(a.data.name, 10), parseInt(b.data.name, 10));
  });
  var svg = container.append('svg:svg')
      .attr('width', w + m[1] + m[3])
      .attr('height', h + m[0] + m[2]);
  svg.append('svg:defs')
    .append('svg:marker')
    .attr('id', 'arrowhead')
    .attr('orient', 'auto')
    .attr('markerWidth', 10)
    .attr('markerHeight', 10)
    .attr('markerUnits', 'strokeWidth')
    .attr('refX', 0)
    .attr('refY', 3)
    .append('svg:path')
    .attr('d', 'M0,0 L0,6 L9,3 z');
  var vis = svg.append('svg:g')
      .attr('transform', 'translate(' + m[3] + ',' + m[0] + ')');

  // Update the nodes.
  var node = vis.selectAll('g.node')
      .data(root.descendants(), function(d) { return d.data.name; });

  var self = this;
  // Enter any new nodes at the parent's previous position.
  var nodeEnter = node.enter().append('svg:g');
  nodeEnter.attr('class', 'population')
    .attr('transform', function(d) { return 'translate(' + d.y + ',' + d.x + ')'; });

  if(NODE_TYPE === 'circle') {
    nodeEnter.append('svg:circle')
      .attr('class', 'outline')
      .attr('r', function(d) { return d.data.radius; });
  } else {
    nodeEnter.append('svg:rect')
      .attr('class', 'outline')
      .attr('x', function(d) { return -d.data.radius; })
      .attr('y', function(d) { return -d.data.radius; })
      .attr('width', function(d) { return 2*d.data.radius; })
      .attr('height', function(d) { return 2*d.data.radius; });
  }

  var choose_colour = function(node_idx, slice_idx, phi) {
    if(samp_colours !== null) {
      let [samp_name, base_colour] = samp_colours[slice_idx];
      let samp_idx = sampnames.indexOf(samp_name);
      let opacity = phi[samp_idx];
      return Util.blend_colours(base_colour, opacity, bg_colour);
    } else if(pop_colours !== null) {
      return pop_colours[node_idx];
    } else {
      return '#428bca';
    }
  };

  var slices = (samp_colours === null) ? 1 : samp_colours.length;
  var add_slice = function(slice_idx, total_slices) {
    if(!(0 <= slice_idx && slice_idx < total_slices)) {
      throw "wrong number of slices: " + slice_idx + ", " + total_slices;
    }
    nodeEnter.append('svg:path')
      .attr('d', function(d) {
        if (NODE_TYPE === 'circle') {
          var deg = 360 / total_slices;
          // Without this, Chrome will draw nothing when there's only a single slice.
          deg = Math.min(deg, 359.9);
          var start = slice_idx * deg;
          var end = (slice_idx + 1)*deg;
          return describeArc(0, 0, d.data.radius, start, end);
        } else {
          var total_width = 2*d.data.radius;
          var slice_width = total_width / total_slices;
          var x = -d.data.radius + (slice_idx * slice_width);
          return describeRect(x, -d.data.radius, slice_width, 2*d.data.radius);
        }
      }).attr('fill', function(d, i) {
        return choose_colour(d.data.name, slice_idx, d.data.phi);
      });
  };
  for(var idx = 0; idx < slices; idx++) {
    add_slice(idx, slices);
  }
  /*nodeEnter.append('svg:path')
    .attr('class', 'divider')
    .attr('stroke', '#aaa')
    .attr('d', function(d) { return 'M 0 -' + d.data.radius + ' V ' + d.data.radius; });*/

  nodeEnter.append('svg:text')
      .attr('font-size', '30')
      .attr('dominant-baseline', 'central')
      .attr('text-anchor', 'middle')
      .text(function(d) { return self._label_node(d.data.name); });

  // Update the links.
  var link = vis.selectAll('path.link')
      .data(root.links(), function(d) { return d.target.data.name; });

  // Enter any new links at the parent's previous position.
  var arrowhead_width = 14;
  link.enter()
    .insert('svg:path', 'g')
    .attr('marker-end', 'url(#arrowhead)')
    .attr('class', 'link')
    .attr('d', d3.linkHorizontal().source(function(d) {
      return [d.source.y + d.source.data.radius, d.source.x];
    }).target(function(d) {
      return [d.target.y - d.source.data.radius - arrowhead_width, d.target.x];
    }));
}

TreePlotter.prototype._generate_tree_struct = function(parents, phi, root_id) {
  var _add_node = function(node_id, struct) {
    struct.name = node_id;
    struct.radius = Math.sqrt(2000 / Math.PI);
    struct.phi = phi[node_id];

    var children = [];
    for(var idx = 0; idx < parents.length; idx++) {
      if(parents[idx] === node_id) {
        children.push(idx + 1);
      }
    }

    struct.children = [];
    children.forEach(function(child_id) {
      var child = {};
      struct.children.push(child);
      _add_node(child_id, child);
    });
  };

  var root = {};
  _add_node(root_id, root);
  return root;
}

TreePlotter.prototype.plot = function(root_id, parents, phi, sampnames, samp_colours, pop_colours, remove_normal, container) {
  container = d3.select(container).append('div');
  var bg_colour = '#ffffff';

  var root = this._generate_tree_struct(parents, phi, root_id);
  this._draw_tree(root, container, sampnames, samp_colours, pop_colours, bg_colour);
  resize_svg(container.selectAll('svg'));
}

