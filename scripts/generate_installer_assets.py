"""
Generate wizard-image.bmp (164x314) and wizard-small.bmp (55x58)
for Inno Setup from existing high-res branding assets.
"""

from pathlib import Path
from PIL import Image

def generate_assets():
    base_dir = Path(__file__).resolve().parent.parent
    branding_dir = base_dir / "frontend" / "assets" / "branding"
    installer_dir = base_dir / "installer"
    installer_dir.mkdir(parents=True, exist_ok=True)

    src_logo = branding_dir / "logo-icon-3d.png"
    if not src_logo.exists():
        src_logo = branding_dir / "favicon.png"

    logo = Image.open(src_logo).convert("RGBA")

    # 1. Wizard Large Image (164 x 314 px) for Inno Setup left banner
    # Dark theme background (#0b0f19)
    bg_large = Image.new("RGB", (164, 314), (11, 15, 25))
    # Resize logo to fit nicely within width 140px
    target_width = 136
    ratio = target_width / float(logo.width)
    target_height = int(logo.height * ratio)
    logo_resized = logo.resize((target_width, target_height), Image.Resampling.LANCZOS)
    
    # Position centered horizontally, slightly towards top
    x_pos = (164 - target_width) // 2
    y_pos = (314 - target_height) // 2 - 20
    bg_large.paste(logo_resized, (x_pos, max(10, y_pos)), logo_resized)
    wizard_large_path = installer_dir / "wizard-image.bmp"
    bg_large.save(wizard_large_path, format="BMP")
    print(f"Generated: {wizard_large_path} ({bg_large.size})")

    # 2. Wizard Small Image (55 x 58 px) for Inno Setup header top-right
    bg_small = Image.new("RGB", (55, 58), (11, 15, 25))
    target_sm_size = 48
    logo_sm = logo.resize((target_sm_size, target_sm_size), Image.Resampling.LANCZOS)
    x_sm = (55 - target_sm_size) // 2
    y_sm = (58 - target_sm_size) // 2
    bg_small.paste(logo_sm, (x_sm, y_sm), logo_sm)
    wizard_small_path = installer_dir / "wizard-small.bmp"
    bg_small.save(wizard_small_path, format="BMP")
    print(f"Generated: {wizard_small_path} ({bg_small.size})")

if __name__ == "__main__":
    generate_assets()
