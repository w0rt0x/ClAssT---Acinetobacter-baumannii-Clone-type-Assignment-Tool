from WebApp import app

if __name__ == '__main__':
    app.run(port=666, debug=False, threaded=True)  # Threading makes Multiple Users at once possible
    # host= '0.0.0.0'

# http://127.0.0.1:666/

# Deployment:
# https://flask.palletsprojects.com/en/master/deploying/
