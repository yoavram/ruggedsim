import os
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

def pandoc():
	print os.system("pandoc -s main.md -o main.pdf --toc --variable=geometry:a4paper")

class PandocEventHandler(FileSystemEventHandler):
	def on_modified(self, event):
		if event.is_directory:
			return
		if event.src_path.endswith('main.md'):
			pandoc()

if __name__ == "__main__":
	pandoc()
	observer = Observer()
	event_handler = PandocEventHandler()
	observer.schedule(event_handler, path='.')
	observer.start()
	try:
		while True:
			time.sleep(1)
	except KeyboardInterrupt:
		observer.stop()
	observer.join()