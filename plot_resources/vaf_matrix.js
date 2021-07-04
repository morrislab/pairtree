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

      var targets = container.querySelectorAll('.matrix tbody tr.' + toggle_type);
      var is_active = E.classList.contains('active');
      if(is_active) {
        E.classList.remove('active');
        E.classList.remove('btn-primary');
        E.classList.add('btn-outline-primary');
        var new_display = 'none';
      } else {
        E.classList.add('active');
        E.classList.remove('btn-outline-primary');
        E.classList.add('btn-primary');
        var new_display = '';
      }
      for(let target of targets) {
        target.style.display = new_display;
      }
    });
  });
}

VafMatrix.prototype._filter_rows = function(container, targets) {
  var rows = container.querySelectorAll('.matrix tbody tr');
  if(targets.length === 0) {
    for(let row of rows) { row.style.display = ''; }
    return;
  }

  rows.forEach((row) => {
    var tid = row.querySelector('.id').textContent.trim();
    if(targets.includes(tid)) {
      row.style.display = '';
    } else {
      row.style.display = 'none';
    }
  });
}

VafMatrix.prototype._configure_filter = function(container) {
  var filter = container.querySelector('.filter');
  var self = this;
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

