"""Apply APA 7 styling to a pandoc-produced docx by editing word/styles.xml in place:
Times New Roman 12 pt throughout, double spacing, black bold headings (L1 centered),
0.5-inch first-line indent for body text, 1-inch margins.
"""
import re
import shutil
import sys
import zipfile
from pathlib import Path

W = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

def restyle(path: Path) -> None:
    tmp = path.with_suffix(".tmp.docx")
    with zipfile.ZipFile(path) as zin:
        names = zin.namelist()
        data = {n: zin.read(n) for n in names}

    styles = data["word/styles.xml"].decode("utf-8")

    # ---- document defaults: TNR 12pt, double spacing --------------------------
    styles = re.sub(
        r"<w:docDefaults>.*?</w:docDefaults>",
        '<w:docDefaults>'
        '<w:rPrDefault><w:rPr>'
        '<w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman" w:eastAsia="Times New Roman" w:cs="Times New Roman"/>'
        '<w:sz w:val="24"/><w:szCs w:val="24"/><w:lang w:val="en-US"/>'
        '</w:rPr></w:rPrDefault>'
        '<w:pPrDefault><w:pPr>'
        '<w:spacing w:after="0" w:line="480" w:lineRule="auto"/>'
        '</w:pPr></w:pPrDefault>'
        '</w:docDefaults>',
        styles, flags=re.S)

    def style_block(sid):
        m = re.search(rf'<w:style [^>]*w:styleId="{sid}"[^>]*>.*?</w:style>', styles, re.S)
        return m

    def replace_style(sid, new):
        nonlocal styles
        m = style_block(sid)
        if m:
            styles = styles[:m.start()] + new + styles[m.end():]

    def heading(sid, name, based, center=False, italic=False, size=24):
        jc = '<w:jc w:val="center"/>' if center else ""
        it = "<w:i/><w:iCs/>" if italic else ""
        return (f'<w:style w:type="paragraph" w:styleId="{sid}">'
                f'<w:name w:val="{name}"/><w:basedOn w:val="{based}"/>'
                f'<w:next w:val="BodyText"/><w:qFormat/>'
                f'<w:pPr><w:keepNext/><w:keepLines/>'
                f'<w:spacing w:before="240" w:after="0" w:line="480" w:lineRule="auto"/>{jc}'
                f'<w:outlineLvl w:val="{int(sid[-1])-1}"/></w:pPr>'
                f'<w:rPr><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/>'
                f'<w:b/><w:bCs/>{it}<w:color w:val="000000"/>'
                f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/></w:rPr></w:style>')

    replace_style("Heading1", heading("Heading1", "heading 1", "Normal", center=True))
    replace_style("Heading2", heading("Heading2", "heading 2", "Normal"))
    replace_style("Heading3", heading("Heading3", "heading 3", "Normal", italic=True))
    replace_style("Heading4", heading("Heading4", "heading 4", "Normal", italic=True))

    # Body Text: first-line indent 0.5" (720 twips), double spaced
    replace_style("BodyText",
        '<w:style w:type="paragraph" w:styleId="BodyText">'
        '<w:name w:val="Body Text"/><w:basedOn w:val="Normal"/><w:qFormat/>'
        '<w:pPr><w:spacing w:after="0" w:line="480" w:lineRule="auto"/>'
        '<w:ind w:firstLine="720"/></w:pPr></w:style>')

    # First Paragraph (pandoc's style after headings): keep APA-consistent indent
    replace_style("FirstParagraph",
        '<w:style w:type="paragraph" w:styleId="FirstParagraph">'
        '<w:name w:val="First Paragraph"/><w:basedOn w:val="BodyText"/><w:qFormat/></w:style>')

    # Compact / table text single-spaced so tables stay readable
    replace_style("Compact",
        '<w:style w:type="paragraph" w:styleId="Compact">'
        '<w:name w:val="Compact"/><w:basedOn w:val="Normal"/>'
        '<w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr></w:style>')

    # Captions: single-spaced, no indent
    replace_style("ImageCaption",
        '<w:style w:type="paragraph" w:styleId="ImageCaption">'
        '<w:name w:val="Image Caption"/><w:basedOn w:val="Normal"/>'
        '<w:pPr><w:spacing w:before="120" w:after="240" w:line="240" w:lineRule="auto"/></w:pPr></w:style>')

    data["word/styles.xml"] = styles.encode("utf-8")

    # ---- 1-inch margins + US Letter in document.xml ---------------------------
    doc = data["word/document.xml"].decode("utf-8")
    doc = re.sub(r'<w:pgSz [^/]*/>',
                 '<w:pgSz w:w="12240" w:h="15840"/>', doc)
    doc = re.sub(r'<w:pgMar [^/]*/>',
                 '<w:pgMar w:top="1440" w:right="1440" w:bottom="1440" w:left="1440" '
                 'w:header="720" w:footer="720" w:gutter="0"/>', doc)
    data["word/document.xml"] = doc.encode("utf-8")

    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for n in names:
            zout.writestr(n, data[n])
    shutil.move(tmp, path)
    print(f"restyled: {path.name}")

for f in sys.argv[1:]:
    restyle(Path(f))
