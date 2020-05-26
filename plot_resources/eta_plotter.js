function EtaPlotter() {
  this._bar_width = 50;
  this._bar_height = 500;
  this._legend_spacing = 60;
  this._font_size = '16px';
  this._legend_font_size = '14px';
  this._label_padding = 10;
  this._legend_splotch_size = 30;
  this._legend_splotch_padding = 5;
  this._legend_splotch_spacing = 20;
  this._col_space = 10;
}

EtaPlotter.prototype._calc_label_width = function(labels) {
  var max_length = 0;
  var char_width = 15;
  for(let label of labels) {
    if(label.length > max_length) {
      max_length = label.length;
    }
  }
  return char_width * max_length;
}

EtaPlotter.prototype._calc_cum = function(mat) {
  var K = mat.length;
  var S = mat[0].length;
  var S_range = Array.from(Array(S).keys());

  var cum = [S_range.map(d => 0)];
  for(let k = 1; k < K; k++) {
    cum.push(S_range.map(s => cum[k - 1][s] + mat[k - 1][s]));
  }
  return cum;
}

EtaPlotter.prototype._plot_etas = function(svg, eta, samp_labels, col_spacing, pop_colours, col_label_height) {
  let self = this;
  let K = eta.length;
  let S = eta[0].length;
  let K_range = Array.from(Array(K).keys());
  let S_range = Array.from(Array(S).keys());
  let eta_cum = this._calc_cum(eta);

  if(samp_labels.length !== S) {
    throw "Wrong number of samp labels";
  }

  let cum_col_spacing = [0].concat(col_spacing);
  for(var idx = 1; idx < cum_col_spacing.length; idx++) {
    cum_col_spacing[idx] += cum_col_spacing[idx - 1];
  }

  let cl = svg.append('svg:g')
    .attr('transform', 'translate(' + (0.5 * self._bar_width) + ',' + (col_label_height - self._label_padding) + ')')
    .selectAll('text')
    .data(samp_labels)
    .join('svg:text')
    .attr('transform', function(d, i) { return 'translate(' + (i*self._bar_width + cum_col_spacing[i]) + ',0) rotate(270)'; })
    .attr('x', 0)
    .attr('y', 0)
    .attr('font-size', this._font_size)
    .attr('font-weight', 'bold')
    .text(function(d, i) { return d; });

  let cols = svg.selectAll('g.col')
    .data(S_range)
    .join('svg:g')
    .attr('class', 'col')
    .attr('transform', function(d, i) {
      return 'translate(' + (i*self._bar_width + cum_col_spacing[i]) + ',' + col_label_height + ')';
    });

  cols.selectAll('rect')
    .data(function(sidx) { return K_range.map(function(k) { return {k: k, s: sidx}; }); })
    .join('svg:rect')
    .attr('width', self._bar_width)
    .attr('height', function(d) { return eta[d.k][d.s] * self._bar_height; })
    .attr('x', 0)
    .attr('y', function(d, i) { return eta_cum[d.k][d.s] * self._bar_height; })
    .attr('fill-opacity', function(d) { return 1.0; })
    .attr('fill', function(d, i) { return pop_colours[i]; });
}

EtaPlotter.prototype._add_pop_legend = function(svg, pop_labels, pop_colours, x_offset, y_offset) {
  let self = this;
  let legend = svg.append('svg:g')
    .attr('transform', 'translate(' + x_offset + ',' + y_offset + ')')
    .selectAll('g.label')
    .data(pop_labels)
    .join('svg:g')
    .attr('transform', function(d, i) { return 'translate(0,' + i*self._legend_splotch_size + ')'; })
    .attr('class', 'label');
  legend
    .append('svg:rect')
    .attr('width', this._legend_splotch_size)
    .attr('height', this._legend_splotch_size)
    .attr('x', 0)
    .attr('y', -0.5*this._legend_splotch_size)
    .attr('fill-opacity', 1.0)
    .attr('fill', function(d, i) { return pop_colours[i]; });
  legend
    .append('svg:text')
    .attr('font-size', this._legend_font_size)
    .attr('x', this._legend_splotch_size + this._legend_splotch_padding)
    .attr('dominant-baseline', 'central')
    .text(function(d) { return d; });
}

EtaPlotter.prototype._remove_small_pops = function(eta, pop_labels, threshold) {
  let K = eta.length;
  let S = eta[0].length;

  for(let k = K - 1; k >= 0; k--) {
    let remove_k = true;
    for(let s = 0; s < S; s++) {
      if(eta[k][s] >= threshold) {
        remove_k = false;
        break;
      }
    }
    if(remove_k) {
      eta.splice(k, 1);
      pop_labels.splice(k, 1);
    }
  }
}

EtaPlotter.prototype.plot = function(eta, samp_labels, container, remove_small_pop_threshold=0.01) {
  let self = this;
  let pop_labels =  Array.from(Array(eta.length).keys()).map(idx => 'Pop. ' + idx);
  if(remove_small_pop_threshold > 0) {
    this._remove_small_pops(eta, pop_labels, remove_small_pop_threshold);
  }

  let K = eta.length;
  let S = eta[0].length;
  let K_range = Array.from(Array(K).keys());
  let S_range = Array.from(Array(S).keys());

  let pop_colours = ColourAssigner.assign_colours(K);
  let pop_label_width = this._calc_label_width(pop_labels);
  let col_label_height = this._calc_label_width(samp_labels);
  let col_spacing = S_range.slice(0, -1).map(s => self._col_space);
  let total_col_spacing = col_spacing.reduce((sum, cur) => sum + cur, 0);

  let legend_x_offset = S*this._bar_width + this._legend_splotch_spacing + total_col_spacing;
  let legend_y_offset = col_label_height + 0.5*this._legend_splotch_size;
  let legend_width = this._legend_splotch_size + this._legend_splotch_padding + pop_label_width;

  let canvas_width = legend_x_offset + legend_width;
  let canvas_height = Math.max(
    col_label_height + this._label_padding + this._bar_height,
    legend_y_offset + K*this._legend_splotch_size,
  );
  let svg = d3.select(container).append('svg:svg')
    .attr('width', canvas_width)
    .attr('height', canvas_height);

  this._plot_etas(
    svg,
    eta,
    samp_labels,
    col_spacing,
    pop_colours,
    col_label_height
  );
  this._add_pop_legend(
    svg,
    pop_labels,
    pop_colours,
    legend_x_offset,
    legend_y_offset
  );
}
