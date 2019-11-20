function PosteriorSumm() {
}

PosteriorSumm.prototype.plot = function(structs, container) {
  if(structs.length === 0) {
    return;
  }

  structs.forEach(function(struct, idx) {
    var row = d3.select(container).append('tr');
    var struct_id = 'structure_' + idx;

    row.append('td').text(idx);
    row.append('td').text(struct.prob.toFixed(3));
    row.append('td').text(struct.nlglh.toFixed(3));
    row.append('td').attr('id', struct_id);

    var root = 0;
    (new TreePlotter()).plot(root, struct.parents, struct.phi, struct.samples, '#' + struct_id);
  });

  resize_svgs();
}

function resize_svgs() {
  document.querySelectorAll('svg').forEach(function(svg) {
    var box = svg.getBBox();
    // Padding is necessary, as `getBBox` doesn't account for borders around
    // nodes -- these would otherwise be clipped off.
    var padding = 6;
    var viewbox = [box.x - 0.5*padding, box.y - 0.5*padding, box.width + padding, box.height + padding];
    svg.setAttribute('viewBox', viewbox.join(' '));
    svg.setAttribute('width', viewbox[2]);
    svg.setAttribute('height', viewbox[3]);
  });
}
