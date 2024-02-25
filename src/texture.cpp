#include <nori/texture.h>
#include <nori/bitmap.h>
NORI_NAMESPACE_BEGIN

class Texture : public Texture_abtract
{
    public:

        Texture(const PropertyList &props)
        {
            loadTexture(props.getString("filename").c_str());
        }
        

        Color3f getPixel(Point2f uv) const
        {
            if(bitmap == nullptr) {
                return defaultVal;
            }
            int imageU = (int)(uv.x() * imageWidht - 0.5) % imageWidht;
            int imageV = (int)(uv.y() * imageHeight - 0.5) % imageHeight;
            
            return bitmap->coeffRef(imageV, imageU);
        }

        void loadTexture(const char* filename)
        {
            Bitmap *textureMap = new Bitmap(filename);
            imageWidht = (int)textureMap->cols();
            imageHeight = (int)textureMap->rows();
            texture = new std::vector<Color3f>();
            texture->reserve(imageWidht * imageHeight);
            for(int y=0;y<imageHeight;y++) {
                for(int x=0;x<imageWidht;x++) {
                    texture->push_back(textureMap->coeffRef(y, x));
                }
            }
            bitmap = textureMap;
        }

        float pdf() const
        {
            return 1.0f / (imageHeight * imageWidht);
        }

        std::string toString() const
        {
            return "Mesh_Texture[]";
        }

    private:
        std::vector<Color3f>* texture;
        int imageWidht;
        int imageHeight;
        Bitmap *bitmap;
        Color3f defaultVal = Color3f(0.0f);
        
};
NORI_REGISTER_CLASS(Texture, "texture");
NORI_NAMESPACE_END