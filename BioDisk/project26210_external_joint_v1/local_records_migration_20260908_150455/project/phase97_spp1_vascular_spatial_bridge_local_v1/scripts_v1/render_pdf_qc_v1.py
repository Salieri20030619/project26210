from pathlib import Path
import pypdfium2 as pdfium
src = Path(r"C:\Users\31937\project26210_spp1_vascular_bridge_v1\figure_v1\PROJECT26210_SPP1_VASCULAR_SPATIAL_BRIDGE_v1.pdf")
out = Path(r"C:\Users\31937\project26210_spp1_vascular_bridge_v1\figure_v1\PROJECT26210_SPP1_VASCULAR_SPATIAL_BRIDGE_v1_pdf_qc.png")
doc = pdfium.PdfDocument(src)
page = doc[0]
bitmap = page.render(scale=2.0)
bitmap.to_pil().save(out)
print(out)
