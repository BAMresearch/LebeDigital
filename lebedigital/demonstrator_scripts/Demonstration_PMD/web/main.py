import subprocess
import os
import http.server
import socketserver
import webbrowser

# start a local fuseki server
def start_fuseki():
    fuseki_path = '' # Update this path
    fuseki_command = os.path.join(fuseki_path, "fuseki-server")
    subprocess.Popen([fuseki_command])
    # upload the datafile ? that was exported by the workflow function

# starting local Sparklis to access the local triple store
def start_webserver():
    PORT = 8000
    Handler = http.server.SimpleHTTPRequestHandler
    with socketserver.TCPServer(("", PORT), Handler) as httpd:
        print(f"Serving at port {PORT}")
        webbrowser.open(f'http://localhost:{PORT}/osparklis.html')  # Open osparklis.html in browser
        httpd.serve_forever()

if __name__ == "__main__":
    start_fuseki()
    start_webserver()
