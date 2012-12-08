#ifndef _PROPERTY_TREE_
#define _PROPERTY_TREE_

#include <string>

#include <json/json.h>

class PropertyTree_jsoncpp {
    public:
        explicit PropertyTree_jsoncpp(const char* conffilename);

        int    asInt()    const { return jsonvalue.asInt();    }
        double asDouble() const { return jsonvalue.asDouble(); }
        std::string asString() const { return jsonvalue.asString(); }

        size_t size() const { return jsonvalue.size(); }

        PropertyTree_jsoncpp operator[] (int i) const {
            return PropertyTree_jsoncpp(jsonvalue[i]);
        }

        bool isMember(const std::string& key) const {
            return jsonvalue.isMember(key);
        }

        PropertyTree_jsoncpp operator[] (const std::string& key) const {
            return PropertyTree_jsoncpp(jsonvalue[key]);
        }

    private:
        explicit PropertyTree_jsoncpp(const Json::Value& value) :
            jsonvalue(value) {}

        Json::Value jsonvalue;
};

typedef PropertyTree_jsoncpp PropertyTree;

#endif // _PROPERTY_TREE_
