// 'letters' or 'digits'
var LABEL_TYPE = 'digits';
// 'circle' or 'square'
var NODE_TYPE = 'circle';

// Following two functions from https://stackoverflow.com/a/18473154.
function polarToCartesian(centerX, centerY, radius, angleInDegrees) {
  var angleInRadians = (angleInDegrees-90) * Math.PI / 180.0;

  return {
    x: centerX + (radius * Math.cos(angleInRadians)),
    y: centerY + (radius * Math.sin(angleInRadians))
  };
}

function describeArc(x, y, radius, startAngle, endAngle) {
  var start = polarToCartesian(x, y, radius, endAngle);
  var end = polarToCartesian(x, y, radius, startAngle);

  var largeArcFlag = endAngle - startAngle <= 180 ? "0" : "1";

  var d = [
    "M", start.x, start.y,
    "A", radius, radius, 0, largeArcFlag, 0, end.x, end.y
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

TreePlotter.colours = {};
TreePlotter.colours.node = '#428bca';
TreePlotter.colours.node_bg = '#ffffff';

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

TreePlotter.prototype._draw_tree = function(root, container, num_pops, left_sample, right_sample) {
  // horiz_padding should be set to the maximum radius of a node, so a node
  // drawn on a boundry won't go over the canvas edge. Since max_area = 8000,
  // we have horiz_padding = sqrt(8000 / pi) =~ 51.
  var horiz_padding = 51;
  var max_depth = this._calculate_max_depth(root);
  var m = [10, horiz_padding, 10, horiz_padding],
      w = Math.max(150, 120*max_depth - m[1] - m[3]),
      h = 430 - m[0] - m[2],
      i = 0;
  var colours = ColourAssigner.assign_colours(num_pops);

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
    .attr('fill', function(d, i) { return colours[d.data.name]; })
    .attr('transform', function(d) { return 'translate(' + d.y + ',' + d.x + ')'; });

  var add_half = function(side) {
    nodeEnter.append('svg:path')
      .attr('class', side + '_half')
      .attr('d', function(d) {
        if (NODE_TYPE === 'circle') {
          if(side === 'left') {
            return describeArc(0, 0, d.data.radius, 180, 360);
          } else {
            return describeArc(0, 0, d.data.radius, 0, 180);
          }
        } else {
          if(side === 'left') {
            return describeRect(-d.data.radius, -d.data.radius, d.data.radius, 2*d.data.radius);
          } else {
            return describeRect(0, -d.data.radius, d.data.radius, 2*d.data.radius);
          }
        }
      }).attr('fill', function(d) {
        var base_colour = TreePlotter.colours.node;
        return Util.rgba2hex(base_colour, d.data.opacities[side], TreePlotter.colours.node_bg);
      });
  };

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
      .attr('height', function(d) { return 2 * d.data.radius; });
  }
  add_half('left');
  add_half('right');
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

TreePlotter.prototype._generate_tree_struct = function(parents, phi, root_id, sampnames, left_sample, right_sample) {
  var left_index = sampnames.indexOf(left_sample);
  var right_index = sampnames.indexOf(right_sample);
  if(left_index === -1 || right_index === -1) {
    throw "Can't find " + left_sample + " or " + right_sample + " samples in " + sampnames;
  }

  var _add_node = function(node_id, struct) {
    struct.name = node_id;
    struct.radius = Math.sqrt(2000 / Math.PI);
    struct.opacities = {
      'left': phi[node_id][left_index],
      'right': phi[node_id][right_index],
    };

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

TreePlotter.prototype.plot = function(root, parents, phi, sampnames, container) {
  var K = phi.length;
  container = d3.select(container).append('div');

  var left_sample = sampnames[0];
  var right_sample = sampnames[0];
  var root = this._generate_tree_struct(parents, phi, root, sampnames, left_sample, right_sample);
  this._draw_tree(root, container, K, left_sample, right_sample);
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
    rgb[channel] = rgb[channel].toString(16).padStart(2, "0");
  });
  var hex = '#' + rgb.r + rgb.g + rgb.b;
  return hex;
}

function PhiMatrix() {
}

PhiMatrix.prototype._calc_ccf = function(phi) {
  var ccf = [];
  // Set clonalidx=0 if you want the phi matrix plots to be unnormalized.
  // Set clonalidx=1 if you want the phi matrix plots to be normalized relative
  // to population 1's phi. (This is only sensible if you don't have any
  // polyclonal trees.)
  var clonalidx = 0;

  for(var rowidx = clonalidx; rowidx < phi.length; rowidx++) {
    ccf[rowidx - 1] = [];
    for(var sampidx = 0; sampidx < phi[rowidx].length; sampidx++) {
      ccf[rowidx - 1].push(phi[rowidx][sampidx] / phi[clonalidx][sampidx]);
    }
  }

  return ccf;
}

PhiMatrix.prototype.plot = function(phi, sampnames, container, convert_to_ccf) {
  if(convert_to_ccf) {
    phi = this._calc_ccf(phi);
  }

  var popnames = phi.map(function(d, i) { return 'Pop. ' + i; });
  var sampcolours = sampnames.map(function(sampname) {
    return '#000000';
  });
  var popcolours = ColourAssigner.assign_colours(phi.length);

  (new MatrixBar()).plot(phi, popnames, popcolours, sampnames, sampcolours, container);
}

function PhiErrorMatrix() {
}

PhiErrorMatrix.prototype.calc_error = function(phi, phi_hat) {
  var error = [];
  for(var i = 0; i < phi.length; i++) {
    error.push([]);
    for(var j = 0; j < phi[i].length; j++) {
      error[i].push(Math.abs(phi[i][j] - phi_hat[i][j]));
    }
  }
  return error;
}

PhiErrorMatrix.prototype._calc_total_error = function(error) {
  var total = 0;
  for(var i = 0; i < error.length; i++) {
    for(var j = 0; j < error[i].length; j++) {
      total += error[i][j];
    }
  }
  return total;
}

PhiErrorMatrix.prototype.plot = function(phi, phi_hat, sampnames, container) {
  var error = this.calc_error(phi, phi_hat);
  d3.select(container).append('h3').text('Total error: ' + this._calc_total_error(error).toFixed(2));
  (new PhiMatrix()).plot(error, sampnames, container);
}

function PhiInterleavedMatrix() {
}

PhiInterleavedMatrix.prototype.plot = function(phi, phi_hat, sampnames, container) {
  var error = (new PhiErrorMatrix()).calc_error(phi, phi_hat);
  var pop_colours = ColourAssigner.assign_colours(phi.length);

  var interleaved = [];
  var row_labels = [];
  var row_colours = [];
  for(var i = 0; i < phi.length; i++) {
    interleaved.push(phi_hat[i]);
    interleaved.push(phi[i]);
    interleaved.push(error[i]);

    row_labels.push('Pop. ' + i + ' data');
    row_labels.push('Pop. ' + i + ' tree');
    row_labels.push('Pop. ' + i + ' error');

    row_colours.push(pop_colours[i]);
    row_colours.push(pop_colours[i]);
    row_colours.push('#000000');
  }

  var sampcolours = sampnames.map(function(sampname) {
    return '#000000';
  });

  (new MatrixBar()).plot(interleaved, row_labels, row_colours, sampnames, sampcolours, container);
}

function MatrixBar() {
}

MatrixBar.prototype._calc_label_width = function(labels) {
  var max_length = 0;
  var char_width = 15;
  for(let label of labels) {
    if(label.length > max_length) {
      max_length = label.length;
    }
  }
  return char_width * max_length;
}

MatrixBar.prototype.plot = function(mat, row_labels, row_colours, col_labels, col_label_colours, container) {
  var num_rows = mat.length;
  var num_cols = mat[0].length;
  var cell_size = 50;
  var row_label_width = this._calc_label_width(row_labels);
  var col_label_height = this._calc_label_width(col_labels);
  var font_size = '24px';
  var label_padding = 10;

  if(row_labels.length !== num_rows) {
    throw "Wrong number of row labels";
  }
  if(col_labels.length !== num_cols) {
    throw "Wrong number of col labels";
  }

  var svg = d3.select(container).append('svg:svg')
    .attr('width', row_label_width + num_cols * cell_size)
    .attr('height', col_label_height + num_rows * cell_size);
  var cl = svg.append('svg:g')
    .attr('transform', function(d, i) { return 'translate(' + (row_label_width + 0.5 * cell_size) + ',' + (col_label_height - label_padding) + ')'; })
    .selectAll('text')
    .data(col_labels)
    .enter()
    .append('svg:text')
    .attr('transform', function(d, i) { return 'translate(' + i * cell_size + ',0) rotate(270)'; })
    .attr('x', 0)
    .attr('y', 0)
    .attr('font-size', font_size)
    .attr('font-weight', 'bold')
    .text(function(d, i) { return d; });
  if(typeof col_label_colours !== 'undefined' && col_label_colours.length === num_cols) {
    cl.attr('fill', function(d, i) { return col_label_colours[i]; });
  }

  var rows = svg.selectAll('g.rows')
    .data(mat)
    .enter()
    .append('svg:g')
    .attr('class', 'rows')
    .attr('transform', function(d, i) { return 'translate(' + row_label_width + ',' + (col_label_height + (i * cell_size)) + ')'; });
  if(typeof row_colours !== 'undefined' && row_colours.length === num_rows) {
    rows.attr('fill', function(d, i) { return row_colours[i]; });
  }

  rows.append('text')
    .attr('x', -label_padding)
    .attr('y', 0.5 * cell_size)
    .attr('dominant-baseline', 'middle')
    .attr('text-anchor', 'end')
    .attr('font-size', font_size)
    .attr('font-weight', 'bold')
    .text(function(d, i) { return row_labels[i]; });
  rows.selectAll('rect')
    .data(function(d) { return d; })
    .enter()
    .append('svg:rect')
    .attr('width', cell_size)
    .attr('height', function(d, i) { return d*cell_size; })
    .attr('x', function(d, i) { return i * cell_size; })
    .attr('y', function(d, i) { return (1 - d)*cell_size; })
    .attr('fill-opacity', function(d) { return 1.0; });
}

function ColourAssigner() {
}

ColourAssigner.assign_colours = function(num_colours) {
  var scale = d3.schemeDark2;

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

function VafMatrix(container) {
  container = document.querySelector(container);
  this._configure_toggles(container);
  this._configure_filter(container);
}

VafMatrix.prototype._configure_toggles = function(container) {
  container.querySelectorAll('.vafmatrix_toggles .btn').forEach((btn) => {
    btn.addEventListener('click', (event) => {
      var E = event.target;
      for(let cls of E.classList) {
        if(cls.startsWith('toggle_')) {
          var toggle_type = cls.replace(/^toggle_/, '');
        }
      }

      var targets = container.querySelectorAll('.matrix tr.' + toggle_type);
      var is_active = E.classList.contains('active');
      if(is_active) {
        E.classList.remove('active');
        var new_display = 'none';
      } else {
        E.classList.add('active');
        var new_display = '';
      }
      for(let target of targets) {
        target.style.display = new_display;
      }
    });
  });
}

VafMatrix.prototype._filter_rows = function(container, targets) {
  var rows = container.querySelectorAll('.matrix tr');
  if(targets.length === 0) {
    for(let row of rows) { row.style.display = ''; }
    return;
  }

  rows.forEach((row) => {
    var tid = row.querySelector('.id').text().trim();
    if(targets.includes(tid)) {
      row.style.display = '';
    } else {
      row.style.display = 'none';
    }
  });
}

VafMatrix.prototype._configure_filter = function(container) {
  var filter = container.querySelector('.filter');
  //var self = this;
  filter.addEventListener('keydown', (E) => {
    if(E.which === 13) {
      E.preventDefault();
      var targets = filter.value.trim().split(',');
      targets = targets.map((T) => { return T.trim(); });
      targets = targets.filter((T) => { return T !== ''; });
      self._filter_rows(container, targets);
    }
  });
}
