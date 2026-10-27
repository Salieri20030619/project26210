import numpy as np
from scipy.stats import rankdata,t,norm

def rz(x):
    x = rankdata(np.asarray(x, float))
    sd = x.std(ddof=1)
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)

def partial(x, y, C):
    rx = x - C @ np.linalg.lstsq(C, x, rcond=None)[0]
    ry = y - C @ np.linalg.lstsq(C, y, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])

def fit(arr):
    n = len(arr)
    A = rz(arr[:, 0])
    X = rz(arr[:, 1])
    Y = rz(arr[:, 2])
    C = np.column_stack([np.ones(n)] + [rz(arr[:, j]) for j in range(3, arr.shape[1])])
    DX = np.column_stack([C, A])
    DY = np.column_stack([C, A, X])
    a = float(np.linalg.lstsq(DX, X, rcond=None)[0][-1])
    yc = np.linalg.lstsq(DY, Y, rcond=None)[0]
    b = float(yc[-1])
    cp = float(yc[-2])
    total = float(np.linalg.lstsq(DX, Y, rcond=None)[0][-1])
    ar = partial(A, X, C)
    br = partial(X, Y, DX)
    k = C.shape[1] - 1

    def pv(r, df):
        return float(2 * t.sf(abs(r) * np.sqrt(df / max(1 - r * r, 1e-15)), df))
    return dict(n=n, a=a, b=b, ab=a * b, cprime=cp, total=total, identity_error=total - cp - a * b, a_rho=ar, a_P=pv(ar, n - k - 2), a_k=k, b_rho=br, b_P=pv(br, n - k - 3), b_k=k + 1, condition_X=float(np.linalg.cond(DX)), condition_Y=float(np.linalg.cond(DY)))

def boot(arr, seed, B=5000):
    rng = np.random.default_rng(seed)
    ii = rng.integers(0, len(arr), (B, len(arr)))
    vr = [rankdata(arr[:, j][ii], axis=1) for j in range(arr.shape[1])]
    vr = [(v - v.mean(1)[:, None]) / v.std(1, ddof=1)[:, None] for v in vr]
    A, X, Y = vr[:3]
    C = [np.ones_like(A)] + vr[3:]

    def reg(y, preds):
        d = np.stack(preds, axis=2)
        g = np.einsum('bni,bnj->bij', d, d)
        rhs = np.einsum('bni,bn->bi', d, y)
        return np.einsum('bij,bj->bi', np.linalg.pinv(g), rhs)
    aa = reg(X, C + [A])[:, -1]
    bb = reg(Y, C + [A, X])[:, -1]
    return (aa * bb, ii)

def RE(est, var):
    w = 1 / np.asarray(var)
    est = np.asarray(est)
    mu = np.average(est, weights=w)
    Q = float(np.sum(w * (est - mu) ** 2))
    den = w.sum() - (w * w).sum() / w.sum()
    tau = max(0, (Q - len(est) + 1) / den)
    rw = 1 / (np.asarray(var) + tau)
    m = float(np.average(est, weights=rw))
    se = float(np.sqrt(1 / rw.sum()))
    return (dict(estimate=m, SE=se, low=m - 1.96 * se, high=m + 1.96 * se, P=float(2 * norm.sf(abs(m / se))), tau2=tau, I2=max(0, (Q - len(est) + 1) / Q) * 100 if Q > 0 else 0), rw / rw.sum())