function EtaPlotter() {
  this._col_width = 50;
  this._col_height = 500;
  this._legend_spacing = 60;
  this._font_size = 16;
  this._legend_font_size = 14;
  this._bar_label_font_size = 14;
  this._label_padding = 10;
  this._legend_splotch_size = 30;
  this._legend_padding = 5;
  this._legend_splotch_spacing = 20;
  this._col_space = 10;
  this._diversity_idx_height = 50;
  this._diversity_idx_spacing = 15;
  this._diversity_idx_legend_width = this._legend_splotch_size;
  this._small_pop_threshold = 0.05;
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

EtaPlotter.prototype._calc_cum_on_axis0 = function(mat) {
  var K = mat.length;
  var S = mat[0].length;
  var S_range = Array.from(Array(S).keys());

  var cum = [S_range.map(d => 0)];
  for(let k = 1; k < K; k++) {
    cum.push(S_range.map(s => cum[k - 1][s] + mat[k - 1][s]));
  }
  return cum;
}

EtaPlotter.prototype._calc_cum = function(A) {
  let cum = [0];
  for(var idx = 1; idx < A.length; idx++) {
    cum.push(cum[idx - 1] + A[idx - 1]);
  }
  return cum;
}

EtaPlotter.prototype._plot_etas = function(svg, eta, pop_labels, samp_labels, col_widths, col_spacing, pop_colours, col_label_height, y_offset) {
  let self = this;
  let K = eta.length;
  let S = eta[0].length;
  let K_range = Array.from(Array(K).keys());
  let S_range = Array.from(Array(S).keys());
  let eta_cum = this._calc_cum_on_axis0(eta);
  let cum_col_spacing = this._calc_cum(col_spacing);
  let cum_col_widths = this._calc_cum(col_widths);

  if(samp_labels.length !== S) {
    throw "Wrong number of samp labels";
  }

  // Setting `pointer-events` to `none` is *critical*. Without this, mousemove
  // and other mouse events will be delivered to the tooltip element, rather
  // than the underlying SVG's rect elements, such that the tooltip will often
  // never display, since its opacity is never changed from 0.
  let tooltip = d3.select('body').append('div')
    .attr('class', 'tooltip bs-tooltip-right')
    .style('pointer-events', 'none')
    .attr('role', 'tooltip')
    .style('opacity', 0);
  let tooltip_inner = tooltip.append('div').attr('class', 'tooltip-inner');

  let _make_percent = V => (V < 0.01 ? '< 1%' : Math.round(100*V) + '%');

  let cl = svg.append('svg:g')
    .attr('class', 'col_labels')
    .attr('transform', 'translate(0,' + (col_label_height - self._label_padding) + ')')
    .selectAll('text')
    .data(samp_labels)
    .join('svg:text')
    .attr('transform', function(d, i) { return 'translate(' + (0.5*col_widths[i] + cum_col_widths[i] + cum_col_spacing[i]) + ',0) rotate(270)'; })
    .attr('x', 0)
    .attr('y', 0)
    .attr('font-size', this._font_size)
    .attr('font-weight', 'bold')
    .text(function(d, i) { return d; });

  let cols = svg.append('svg:g')
    .attr('class', 'cols')
    .attr('transform', 'translate(0,' + y_offset + ')')
    .selectAll('g.col')
    .data(S_range)
    .join('svg:g')
    .attr('class', 'col')
    .attr('transform', function(d, i) {
      return 'translate(' + (cum_col_widths[i] + cum_col_spacing[i]) + ',0)';
    });

  cols.selectAll('rect')
    .data(function(sidx) { return K_range.map(function(k) { return {k: k, s: sidx}; }); })
    .join('svg:rect')
    .attr('width', d => col_widths[d.s])
    .attr('height', d => eta[d.k][d.s] * self._col_height)
    .attr('x', 0)
    .attr('y', function(d, i) { return eta_cum[d.k][d.s] * self._col_height; })
    .attr('fill-opacity', function(d) { return 1.0; })
    .attr('fill', function(d, i) { return pop_colours[i]; })
    .on('mousemove', function(d, i) {
      tooltip.style('opacity', 0.9).style('top', (d3.event.pageY - 10) + 'px').style('left', (d3.event.pageX + 10) + 'px');
      tooltip_inner.text(pop_labels[i] + ': ' + _make_percent(eta[d.k][d.s]));
    }).on('mouseout', function(d) {
      tooltip.style('opacity', 0);
    });

  let bar_labels = cols.selectAll('text')
    .data(function(sidx) { return K_range.map(function(k) { return {k: k, s: sidx}; }); })
    .join('svg:text')
    // I subtract 0.5*self._bar_label_font_size here to get the vertical
    // alignment of the two text lines approximately correct. I'm not sure the
    // math is right, but it produces a (mostly?) correct visual result.
    .style('pointer-events', 'none') // Pass events through to underlying `rect`.
    .attr('y', function(d) { return (0.5*eta[d.k][d.s] + eta_cum[d.k][d.s])*self._col_height - 0.5*self._bar_label_font_size; })
    .attr('display', function(d) { return eta[d.k][d.s]*self._col_height >= 2.0*self._bar_label_font_size ? 'inline' : 'none'})
    .attr('font-size', this._bar_label_font_size)
    .attr('dominant-baseline', 'central')
    .attr('text-anchor', 'middle');

  bar_labels.append('tspan')
    .attr('x', d => 0.5*col_widths[d.s])
    .text(function(d, i) { return pop_labels[i] });
  bar_labels.append('tspan')
    .attr('x', d => 0.5*col_widths[d.s])
    .attr('dy', 1.2*this._bar_label_font_size)
    .text((d, i) => _make_percent(eta[d.k][d.s]));
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
    .attr('x', this._legend_splotch_size + this._legend_padding)
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

EtaPlotter.prototype._plot_diversity_indices = function(svg, cdi, cmdi, y_offset, col_widths, col_spacing, colour) {
  let self = this;
  let cum_col_width = this._calc_cum(col_widths);
  let cum_col_spacing = this._calc_cum(col_spacing);

  let container = svg.append('svg:g')
    .attr('class', 'diversity_idx_bars')
    .attr('transform', 'translate(0,' + y_offset + ')');

  let max_cdi = Math.max(...cdi);
  let max_cmdi = Math.max(...cmdi);
  let div_idxs = new Map([
    ['cdi', cdi.map(C => C / max_cdi)],
    ['cmdi', cmdi.map(C => C / max_cmdi)],
  ]);
  let y_offsets = new Map([['cdi', 0], ['cmdi', self._diversity_idx_height + self._diversity_idx_spacing]]);

  div_idxs.forEach((val, key) => {
    container.selectAll('.' + key)
      .data(val)
      .join('svg:rect')
      .attr('class', key)
      .attr('x', (d, i) => cum_col_width[i] + cum_col_spacing[i])
      .attr('y', V => (1 - V)*self._diversity_idx_height + y_offsets.get(key))
      .attr('width', (d, i) => col_widths[i])
      .attr('height', V => V*self._diversity_idx_height)
      .attr('fill', V => colour(V));
  });
}

EtaPlotter.prototype._ramp = function(colour, n = 256) {
  const canvas = d3.create('canvas') .attr('width', 1) .attr('height', n)
    .node();
  const context = canvas.getContext("2d");
  for (let i = 0; i < n; ++i) {
    context.fillStyle = colour(1 - i/(n - 1));
    context.fillRect(0, i, 1, 1);
  }
  return canvas;
}

EtaPlotter.prototype._make_diversity_idx_legend = function(svg, colour, x_offset, y_offset) {
  let scale = d3.scaleSequential([0, 100], colour);

  let legend = svg.append('svg:g')
    .attr('class', 'diversity_idx_legend')
    .attr('transform', 'translate(' + x_offset + ',' + y_offset + ')');

  let legend_offsets = new Map([
    ['cdi', 0],
    ['cmdi', this._diversity_idx_height + this._diversity_idx_spacing]
  ]);
  let labels = new Map([
    ['cdi', 'Clone'],
    ['cmdi', 'Clone and mutation'],
  ]);

  legend_offsets.forEach((Y, key) => {
    legend.append('image')
      .attr('x', 0)
      .attr('y', Y)
      .attr('width', this._diversity_idx_legend_width)
      .attr('height', this._diversity_idx_height)
      .attr('preserveAspectRatio', 'none')
      .attr('xlink:href', this._ramp(scale.interpolator()).toDataURL());

    let axis_scale = d3.scaleLinear()
      .domain([0, 1])
      .range([this._diversity_idx_height, 0]);
    let font_size = 0.9*this._legend_font_size;

    let label = legend.append('svg:g');
    label.append('svg:text').text(labels.get(key));
    label.append('svg:text').attr('dy', font_size).text('diversity index');
    label.selectAll('text')
      .attr('x', 0)
      .attr('y', 0)
      .attr('dominant-baseline', 'hanging')
      .attr('font-size', font_size);
    let label_offset = 0.5*(this._diversity_idx_height - label.node().getBoundingClientRect().height);
    label.attr('transform', 'translate(' + (this._diversity_idx_legend_width + 2*this._legend_padding) + ',' + (Y + label_offset) + ')');

    legend.append('g')
      .attr('transform', 'translate(' + (this._diversity_idx_legend_width + this._legend_padding) + ',' + Y + ')')
      .call(d3.axisRight(axis_scale)
        .ticks(1)
        .tickSize(5)
      ).call(g => g.selectAll('text').attr('font-size', font_size));
  });
}

EtaPlotter.prototype._renormalize_eta = function(eta) {
  let K = eta.length;
  let S = eta[0].length;

  for(let s = 0; s < S; s++) {
    let eta_s_sum = eta.reduce((sum, cur) => sum + cur[s], 0);
    for(let k = 0; k < K; k++) {
      eta[k][s] /= eta_s_sum;
    }
  }
}

EtaPlotter.prototype._add_yaxis = function(svg, height, yoffset) {
  let scale = d3.scaleLinear()
    .domain([0, 1])
    .range([height, 0]);
  let axis = d3.axisLeft(scale)
    .tickFormat(x => Math.round(100*x) + '%');
  let container = svg.append('g')
    .attr('transform', 'translate(-10,' + yoffset + ')');
  container.append('text')
    .attr('text-anchor', 'middle')
    .attr('transform', 'rotate(-90)')
    .attr('x', -0.5*height)
    .attr('y', -60)
    .attr('font-size', this._bar_label_font_size)
    .text('Population frequency');
  container.append('g')
    .call(axis)
    .selectAll('text')
    .attr('font-size', this._bar_label_font_size);
}

EtaPlotter.prototype.plot = function(eta, cdi, cmdi, samp_labels, container, remove_small_pops=false, remove_normal=false) {
  let self = this;
  let pop_labels =  Array.from(Array(eta.length).keys()).map(idx => 'Pop. ' + idx);

  if(remove_normal) {
    eta = eta.slice(1);
    pop_labels = pop_labels.slice(1);
  }
  if(remove_small_pops) {
    this._remove_small_pops(eta, pop_labels, this._small_pop_threshold);
  }
  this._renormalize_eta(eta);

  let K = eta.length;
  let S = eta[0].length;
  let K_range = Array.from(Array(K).keys());
  let S_range = Array.from(Array(S).keys());

  this._col_height = Math.max(this._col_height, K*this._legend_splotch_size);
  let pop_colours = ColourAssigner.assign_colours(K);
  let pop_label_width = this._calc_label_width(pop_labels);
  let col_label_height = this._calc_label_width(samp_labels);

  let col_widths = samp_labels.map(label => (label.indexOf('Xeno') > -1 ? self._col_width : 2*self._col_width));
  let total_col_widths = col_widths.reduce((sum, cur) => sum + cur, 0);

  let col_spacing = S_range.map(s => self._col_space);
  col_spacing[S - 1] = 0;
  let total_col_spacing = col_spacing.reduce((sum, cur) => sum + cur, 0);

  let bar_y_offset = 2*this._diversity_idx_height + 2*this._diversity_idx_spacing + col_label_height;
  let legend_x_offset = total_col_widths + this._legend_splotch_spacing + total_col_spacing;
  let legend_y_offset = bar_y_offset + 0.5*this._legend_splotch_size;
  let legend_width = this._legend_splotch_size + this._legend_padding + pop_label_width;

  let svg = d3.select(container).append('svg:svg');

  let diversity_idx_colour = d3.interpolateViridis;
  this._plot_etas(
    svg,
    eta,
    pop_labels,
    samp_labels,
    col_widths,
    col_spacing,
    pop_colours,
    col_label_height,
    bar_y_offset
  );
  this._plot_diversity_indices(
    svg,
    cdi,
    cmdi,
    col_label_height,
    col_widths,
    col_spacing,
    diversity_idx_colour,
  );
  this._add_pop_legend(
    svg,
    pop_labels,
    pop_colours,
    legend_x_offset,
    legend_y_offset
  );
  this._make_diversity_idx_legend(
    svg,
    diversity_idx_colour,
    legend_x_offset,
    col_label_height,
  );
  this._add_yaxis(
    svg,
    this._col_height,
    bar_y_offset
  );

  resize_svg(svg);
}
