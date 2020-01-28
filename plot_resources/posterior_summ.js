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
    row.append('td').text((100 * struct.prob).toFixed(1) + '%');
    row.append('td').text(struct.nlglh.toFixed(3));
    row.append('td').text(struct.count);
    row.append('td').attr('id', struct_id);

    var root = 0;
    (new TreePlotter()).plot(root, struct.parents, struct.phi, struct.samples, '#' + struct_id);
  });
}
