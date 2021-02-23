# coding: utf-8
import numpy as np
from PIL import Image
import astra
from tomopy.misc.phantom import shepp2d

def fp(pg, vg, v):
    proj_id = astra.create_projector('line_fanflat', pg, vg);
    v_id = astra.data2d.create('-vol', vg, v)
    rt_id = astra.data2d.create('-sino', pg)
    fp_cfg = astra.astra_dict('FP')
    fp_cfg['ProjectorId'] = proj_id
    fp_cfg['VolumeDataId'] = v_id
    fp_cfg['ProjectionDataId'] = rt_id
    fp_id = astra.algorithm.create(fp_cfg)
    astra.algorithm.run(fp_id)
    out = astra.data2d.get(rt_id)
    astra.algorithm.delete(fp_id)
    astra.data2d.delete(rt_id)
    astra.data2d.delete(v_id)
    return out

def main():
    size = 256
    angles = np.linspace(0, 2*np.pi, 360, False)
    source_object = 7000
    object_det = 500
    phantom = np.squeeze(shepp2d(size))
    vol_geom = astra.create_vol_geom((size, size))
    proj_geom = astra.create_proj_geom("fanflat", 1.0, size, angles, source_object, object_det)
    sinogram = fp(proj_geom, vol_geom, phantom)
    with open("sinogram.txt", "w") as f :
        for i in range(len(angles)):
            f.write(" ".join([str(x) for x in sigogram[i, :].tolist()]))
            f.write("\n")
    sinogram = sinogram/sinogram.max()*255.0
    im = Image.fromarray(np.transpose(phantom.astype(np.uint8)))
    im.convert("L")
    im.save("phantom_sq.bmp")
    im = Image.fromarray(np.transpose(sinogram.astype(np.uint8)))
    im.convert("L")
    im.save("sinogram_sq.bmp")

main()
