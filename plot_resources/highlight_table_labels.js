function highlight_labels(cell, should_highlight) {
  var row = cell.closest('tr');
  var tbl = cell.closest('table');
  var col_idx = Array.from(row.children).indexOf(cell);
  var col_label = tbl.querySelector('thead tr').children[col_idx];
  var row_label = row.children[0];

  for(let elem of [col_label, row_label]) {
    var cls = 'highlighted';
    if(elem.classList.contains(cls)) {
      elem.classList.remove(cls);
    } else {
      elem.classList.add(cls);
    }
  }
}

document.addEventListener('DOMContentLoaded', () => {
  document.querySelectorAll('table.matrix td').forEach((cell) => {
    cell.addEventListener('mouseenter', (E) => {
      highlight_labels(E.target, true);
      E.preventDefault();
    });
    cell.addEventListener('mouseleave', (E) => {
      highlight_labels(E.target, false);
      E.preventDefault();
    });
  });
});
