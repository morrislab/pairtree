var IS_XENO = window.location.toString().indexOf('pairwise.xeno') > -1;
IS_XENO = true;

function TreePlotter() {
}

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
TreePlotter.colours.diag = '#428bca';
TreePlotter.colours.relapse = '#ca4242';
TreePlotter.colours.other = '#d4b831';
TreePlotter.colours.node_bg = '#ffffff';

TreePlotter.prototype._label_node = function(node_id) {
  // If we're examining a patient+xeno tree, label the node with the node_id.
  // Otherwise, we're presumably examining a patient-only tree, so label it
  // with a letter.
  if(IS_XENO) {
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
        if (IS_XENO) {
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
        // Use other_colour for xenos.
        if(side === 'left') {
          var base_colour = (left_sample.indexOf('X') > -1)  ? TreePlotter.colours.other : TreePlotter.colours.diag;
        } else {
          var base_colour = (right_sample.indexOf('X') > -1) ? TreePlotter.colours.other : TreePlotter.colours.relapse;
        }
        return Util.rgba2hex(base_colour, d.data.opacities[side], TreePlotter.colours.node_bg);
      });
  };

  if(IS_XENO) {
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
      // Uncomment this if you want to disallow polyclonal trees.
      //throw "Unexpected children from root " + summary.trees[tidx].root + ": " + struct[root];
    }
    var clonal = struct[root][0];

    // Note to self: change "root" to "clonal" in the call below to eliminate
    // the normal node 0 when drawing trees.
    var root = self._generate_tree_struct(summary.params.samples, struct, pops, root, left_sample, right_sample);
    self._draw_tree(root, container, Object.keys(pops).length, left_sample, right_sample);
  });
}

function TreeUtil() {
}

TreeUtil.calc_radius = function(scale) {
  var min_area = 700, max_area = 8000;
  // scale must be in [0, 1].
  //var area = min_area + scale*(max_area - min_area);
  var area = 2000;
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

PhiMatrix.prototype._filter_samples = function(phi_matrix, sampnames, samps_to_keep) {
  var filtered_samps = [];
  var filtered_phi = [];
  phi_matrix.forEach(function(row) {
    filtered_phi.push([]);
  });

  samps_to_keep.forEach(function(samp) {
    var sidx = sampnames.indexOf(samp);
    if(sidx < 0) {
      throw "Can't find " + samp + " in " + sampnames;
    }
    filtered_samps.push(samp);
    phi_matrix.forEach(function(phi, pidx) {
      filtered_phi[pidx].push(phi_matrix[pidx][sidx]);
    });
  });

  return {'phi': filtered_phi, 'samples': filtered_samps};
}

PhiMatrix.prototype.plot = function(sampid, phi_path, container) {
  var self = this;
  var filters = {
    //'SJBALL022609': ['D', 'dPDX 26', 'dPDX 8', 'dPDX 2', 'dPDX 15', 'dPDX 7', 'dPDX 14', 'dPDX 20', 'dPDX 29', 'R1', 'rPDX 25', 'rPDX 2', 'rPDX 17', 'rPDX 21']
  };
  d3.json(phi_path, function(phi_data) {
    if(IS_XENO && filters.hasOwnProperty(sampid)) {
      phi_data = self._filter_samples(
        phi_data['phi'],
        phi_data['samples'],
        filters[sampid]
      );
    }

    var ccf = self._calc_ccf(phi_data['phi']);
    var popnames = ccf.map(function(d, i) { return 'Pop. ' + (i + 1); });
    var sampnames = phi_data['samples'];
    var sampcolours = sampnames.map(function(sampname) {
      return sampname[0].toUpperCase() === 'D' ? TreePlotter.colours.diag : TreePlotter.colours.relapse;
    });
    var popcolours = ColourAssigner.assign_colours(ccf.length);

    (new MatrixBar()).plot(ccf, popnames, popcolours, sampnames, sampcolours, container);
  });
}

function MatrixBar() {
}

MatrixBar.prototype.plot = function(mat, row_labels, row_colours, col_labels, col_label_colours, container) {
  var num_rows = mat.length;
  var num_cols = mat[0].length;
  var cell_size = 50;
  var row_label_width = 100;
  var col_label_height = 120;
  var font_size = '24px';
  var label_padding = 10;

  if(row_labels.length !== num_rows) {
    throw "Wrong number of row labels";
  }
  if(col_labels.length !== num_cols) {
    throw "Wrong number of col labels";
  }

  var svg = d3.select(container).html('').append('svg:svg')
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

function ClusterMatrix() {
}

ClusterMatrix.prototype.plot = function(clustermat_path, container) {
  d3.json(clustermat_path, function(J) {
    var clustermat = J.clustermat;
    var row_labels = clustermat.map(function(row, idx) {
      return idx === 0 ? 'Garbage' : 'Pop. ' + idx;
    });
    var col_labels = clustermat[0].map(function(col, idx) {
      var letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
      var popname = letters[idx - 1];
      return idx === 0 ? 'Garbage' : 'Pop. ' + popname;
    });
    var row_colours = ColourAssigner.assign_colours(clustermat.length);
    var col_label_colours = undefined;
    (new MatrixBar()).plot(clustermat, row_labels, row_colours, col_labels, col_label_colours, container);
  });
}

function VafMatrix(container) {
  this._configure_toggles(container);
}

VafMatrix.prototype._configure_toggles = function(container) {
  var togglers = {
    phi: function() {
      return $(this).find('td.id').text().startsWith('P');
    },
    cluster_means: function() {
      return $(this).find('td.id').text().startsWith('C');
    },
    cluster_members: function() {
      return !togglers.garbage.call(this) && $(this).find('td.id').text().startsWith('s');
    },
    garbage: function() {
      var cluster = parseInt($(this).find('td.cluster').text(), 10);
      return isNaN(cluster);
    }
  };

  $(container).find('.vafmatrix_toggles .btn').change(function() {
    var E = $(this);
    var active = E.hasClass('active');
    var toggle_type = E.attr('class').split(/\s+/).filter(function(cls) { return cls.startsWith('toggle_'); })[0].replace(/^toggle_/, '');
    var targets = $(container).find('.matrix tr').filter(togglers[toggle_type]);

    if(active) {
      E.removeClass('active');
      targets.hide();
    } else {
      E.addClass('active');
      targets.show();
    }
  });
}
