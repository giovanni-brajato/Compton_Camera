from PIL import Image

def newImg():
    img = Image.new('RGB', (100, 100))
    img.putpixel((30,60), (155,155,55))
    return img

wallpaper = newImg()
wallpaper.show()