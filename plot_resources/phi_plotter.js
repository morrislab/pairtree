function PhiMatrix() {
}

PhiMatrix.prototype.plot = function(phi, parents, sampnames, container, orientation, remove_normal=false) {
  // Note that if `remove_normal = false`, you can set `parents = null` when
  // calling this function.

  var popnames = phi.map(function(d, i) { return 'Pop. ' + i; });
  var popcolours = ColourAssigner.assign_colours(phi.length);
  var sampcolours = null;

  // By doing this operation after initializing `popnames` and `popcolours`, we
  // keep colours and names consistent regardless of whether `remove_normal` is
  // set.
  if (remove_normal) {
    phi = Util.calc_ccf(phi, parents);
    popnames = popnames.slice(1);
    popcolours = popcolours.slice(1);
  }

  if(orientation === 'samples_as_rows') {
    let phi_T = Util.transpose(phi);
    (new MatrixBar()).plot(phi_T, sampnames, sampcolours, popnames, popcolours, container);
  } else {
    (new MatrixBar()).plot(phi, popnames, popcolours, sampnames, sampcolours, container);
  }
}

function ErrorCalculator() {
}

ErrorCalculator.calc_error = function(phi, phi_hat) {
  var error = [];
  for(var i = 0; i < phi.length; i++) {
    error.push([]);
    for(var j = 0; j < phi[i].length; j++) {
      error[i].push(Math.abs(phi[i][j] - phi_hat[i][j]));
    }
  }
  return error;
}

ErrorCalculator.calc_total_error = function(error) {
  var total = 0;
  for(var i = 0; i < error.length; i++) {
    for(var j = 0; j < error[i].length; j++) {
      total += error[i][j];
    }
  }
  return total;
}

function PhiInterleavedMatrix() {
}

PhiInterleavedMatrix.prototype.plot = function(phi, phi_hat, sampnames, container, orientation) {
  var error = ErrorCalculator.calc_error(phi, phi_hat);
  var total_error = ErrorCalculator.calc_total_error(error);
  d3.select(container).append('h3').text('Total error: ' + total_error.toFixed(2));

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

  var sampcolours = null;
  if(orientation === 'samples_as_rows') {
    let interleaved_T = Util.transpose(interleaved);
    (new MatrixBar()).plot(interleaved_T, sampnames, sampcolours, row_labels, row_colours, container);
  } else {
    (new MatrixBar()).plot(interleaved, row_labels, row_colours, sampnames, sampcolours, container);
  }
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

MatrixBar.prototype.plot = function(mat, row_labels, row_colours, col_labels, col_colours, container) {
  var num_rows = mat.length;
  var num_cols = mat[0].length;
  var cell_size = 50;
  var row_label_width = this._calc_label_width(row_labels);
  var col_label_height = this._calc_label_width(col_labels);
  var font_size = '16px';
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
  if(col_colours && col_colours.length === num_cols) {
    cl.attr('fill', function(d, i) { return col_colours[i]; });
  } else {
    cl.attr('fill', '#000');
  }

  var rows = svg.selectAll('g.rows')
    .data(mat)
    .enter()
    .append('svg:g')
    .attr('class', 'rows')
    .attr('transform', function(d, i) { return 'translate(' + row_label_width + ',' + (col_label_height + (i * cell_size)) + ')'; });

  let rl = rows.append('text')
    .attr('x', -label_padding)
    .attr('y', 0.5 * cell_size)
    .attr('dominant-baseline', 'middle')
    .attr('text-anchor', 'end')
    .attr('font-size', font_size)
    .attr('font-weight', 'bold')
    .text(function(d, i) { return row_labels[i]; });
  if(row_colours && row_colours.length === num_rows) {
    rl.attr('fill', function(d, i) { return row_colours[i]; });
  }

  rows.selectAll('rect')
    .data((d, i) => d.map(function(val, col_idx) { return { val: val, col_idx: col_idx, row_idx: i}; }))
    .enter()
    .append('svg:rect')
    .attr('width', cell_size)
    .attr('height', function(d) { return d.val*cell_size; })
    .attr('x', function(d) { return d.col_idx * cell_size; })
    .attr('y', function(d) { return (1 - d.val)*cell_size; })
    .attr('fill', function(d) {
      // col_colours overrides row_colours.
      if(col_colours && col_colours.length == num_cols) {
        return col_colours[d.col_idx];
      } else if(row_colours && row_colours.length == num_rows) {
        return row_colours[d.row_idx];
      }
      return '#000';
    }).attr('fill-opacity', function(d) { return 1.0; });
}
