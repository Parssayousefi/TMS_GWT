# -*- coding: utf-8 -*-
"""Wrap a PNG in a PDF losslessly. PNG IDAT is deflate-compressed with the
PNG row predictor, which is exactly what PDF's FlateDecode + /Predictor 15
expects, so the pixel data is copied through with NO recompression and NO
quality loss. Notebook code is untouched."""
import struct, zlib, sys
from pathlib import Path
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

def png_chunks(b):
    assert b[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    i = 8
    while i < len(b):
        ln = struct.unpack(">I", b[i:i + 4])[0]
        typ = b[i + 4:i + 8]
        yield typ, b[i + 8:i + 8 + ln]
        i += 12 + ln

def png_to_pdf(png_path, pdf_path, page_w_pt, page_h_pt):
    raw = Path(png_path).read_bytes()
    hdr = idat = None
    for typ, data in png_chunks(raw):
        if typ == b"IHDR": hdr = data
        elif typ == b"IDAT": idat = (idat or b"") + data
    w, h, bits, ctype, comp, filt, inter = struct.unpack(">IIBBBBB", hdr)
    assert comp == 0 and filt == 0 and inter == 0, "unsupported PNG variant"
    assert bits == 8 and ctype in (2, 0), f"need 8-bit RGB/grey, got ctype={ctype}"
    colors = 3 if ctype == 2 else 1
    cs = "/DeviceRGB" if colors == 3 else "/DeviceGray"

    objs = []
    objs.append(b"<< /Type /Catalog /Pages 2 0 R >>")
    objs.append(b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>")
    objs.append((f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {page_w_pt:.4f} "
                 f"{page_h_pt:.4f}] /Resources << /XObject << /Im0 4 0 R >> >> "
                 f"/Contents 5 0 R >>").encode())
    img = (f"<< /Type /XObject /Subtype /Image /Width {w} /Height {h} "
           f"/ColorSpace {cs} /BitsPerComponent 8 /Filter /FlateDecode "
           f"/DecodeParms << /Predictor 15 /Colors {colors} "
           f"/BitsPerComponent 8 /Columns {w} >> /Length {len(idat)} >>").encode()
    objs.append(img + b"\nstream\n" + idat + b"\nendstream")
    content = f"q {page_w_pt:.4f} 0 0 {page_h_pt:.4f} 0 0 cm /Im0 Do Q".encode()
    objs.append(f"<< /Length {len(content)} >>".encode() +
                b"\nstream\n" + content + b"\nendstream")

    out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
    offs = []
    for n, body in enumerate(objs, 1):
        offs.append(len(out))
        out += f"{n} 0 obj\n".encode() + body + b"\nendobj\n"
    xref = len(out)
    out += f"xref\n0 {len(objs)+1}\n".encode() + b"0000000000 65535 f \n"
    for o in offs:
        out += f"{o:010d} 00000 n \n".encode()
    out += (f"trailer\n<< /Size {len(objs)+1} /Root 1 0 R >>\nstartxref\n"
            f"{xref}\n%%EOF\n").encode()
    Path(pdf_path).write_bytes(out)
    return w, h, len(out)

if __name__ == "__main__":
    SP = Path(r"C:\Users\YousefiP\AppData\Local\Temp\claude\D--projects-GWT\0d0eb006-154d-436e-9ee6-c43fba899665\scratchpad")
    F = Path(r"D:\projects\TMS_GWT\pipeline\manuscript\figures_final")
    # flatten alpha onto white and re-save as plain 8-bit RGB PNG
    im = Image.open(F / "FigS7.png")
    if im.mode != "RGB":
        bg = Image.new("RGB", im.size, "white")
        bg.paste(im, mask=im.split()[-1] if im.mode in ("RGBA", "LA") else None)
        im = bg
    flat = SP / "_s7_flat.png"
    im.save(flat, "PNG", optimize=True, compress_level=9)
    print(f"flattened PNG: {flat.stat().st_size/1e6:.2f} MB  {im.size}")
    w, h, n = png_to_pdf(flat, SP / "_s7_lossless.pdf", 1396.8, 1694.030625)
    print(f"lossless PDF : {n/1e6:.2f} MB  ({w}x{h} px)")
