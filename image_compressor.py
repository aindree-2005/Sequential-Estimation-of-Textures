from PIL import Image
img = Image.open("patt.png")
img_50 = img.resize((50, 50), Image.LANCZOS)
img_50.save("pattern_cell.png")
