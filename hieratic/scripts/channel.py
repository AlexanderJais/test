import cv2, numpy as np, math

def _rot(img, ang):
    h,w=img.shape[:2]; d=int(math.hypot(h,w))+4
    M=cv2.getRotationMatrix2D((w/2,h/2),ang,1.0)
    M[0,2]+=d/2-w/2; M[1,2]+=d/2-h/2
    return cv2.warpAffine(img,M,(d,d),flags=cv2.INTER_NEAREST),M

def _unrot(img,M,shape):
    return cv2.warpAffine(img,cv2.invertAffineTransform(M),(shape[1],shape[0]),flags=cv2.INTER_NEAREST)

def scan_gaps(SD,HL,ang):
    A,_=_rot(SD.astype(np.uint8),ang); B,_=_rot(HL.astype(np.uint8),ang)
    W=A.shape[1]; ar=np.arange(W)
    idx=np.where(A>0,ar[None,:],-1)
    last=np.maximum.accumulate(idx,axis=1)
    nxt=np.minimum.accumulate(np.where(B>0,ar[None,:],W+1)[:,::-1],axis=1)[:,::-1]
    m=(A>0)&(nxt<=W)
    return (nxt-ar[None,:])[m]

def channel_between(SD,HL,ang,Wmax):
    """Pixels with a shadow rim behind and a highlight rim ahead along +u,
    closer together than Wmax: the incised channel itself, not its rims."""
    A,M=_rot(SD.astype(np.uint8),ang); B,_=_rot(HL.astype(np.uint8),ang)
    H,W=A.shape; ar=np.arange(W)
    idx=np.where(A>0,ar[None,:],-1)
    last=np.maximum.accumulate(idx,axis=1)                 # nearest shadow at or before x
    back=np.where(last>=0,ar[None,:]-last,10**6)
    nxt=np.minimum.accumulate(np.where(B>0,ar[None,:],10**6)[:,::-1],axis=1)[:,::-1]
    fwd=np.where(nxt<10**6,nxt-ar[None,:],10**6)           # nearest highlight at or after x
    ch=(((back+fwd)<=Wmax)&(back<10**6)&(fwd<10**6)).astype(np.uint8)
    return _unrot(ch,M,SD.shape)
