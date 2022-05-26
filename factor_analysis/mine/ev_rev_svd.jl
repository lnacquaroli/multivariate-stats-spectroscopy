"""Perform thin-SVD."""
function ev_rev_svd(X, rowsX, colsX, sm)
    F = svd(X, full=false) # F.U, F.S, F.V
    ev = (F.S).^2
    rev = [ ev[j]/((rowsX - j + 1 )*(colsX - j + 1)) for j=1:sm ]
    return F, ev, rev
end
