#include <iostream>
#include "httplib.h"
#include "src/Simplify.h"
#include "src/Obj.h"

using namespace httplib;
using namespace std;

int main() {
    int port = 8080;
    Server svr;

    svr.Post("/simplify", [](const Request &req, Response &res) {
        Obj obj;
        cout << "Received a request." << endl;

        string config;
        if (req.has_file("config")) {
            auto file = req.get_file_value("config");
            config = file.content;
            cout << "Receive config string: " << config << endl;
        }

        if (req.has_file("file")) {
            auto file = req.get_file_value("file");

            cout << "File name: " << file.filename << endl;
            cout << "File size: " << file.content.size() << " bytes" << endl;

            //sim.load_obj(file.content, false);

            obj.load(file.content);
            if (config.find("mesh"))
            {
                Simplify sim(obj);
                sim.load_obj();
                sim.simplify_mesh(0.04, 7.0, true);
                cout << "Simplify success!" << endl;
            }

            string model_data = obj.write_to_string();
            res.status = 200;
            res.set_content(model_data, "application/octet-stream");
        }
    });
    svr.listen("0.0.0.0", port);
    return 0;
}