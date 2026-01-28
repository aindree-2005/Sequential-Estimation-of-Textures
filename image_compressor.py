from PIL import Image
img = Image.open("weave.png")
img_50 = img.resize((100,100), Image.LANCZOS)
img_50.save("weave_comp.png")