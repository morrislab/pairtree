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
