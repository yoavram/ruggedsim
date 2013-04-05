import os
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

class PandocEventHandler(FileSystemEventHandler):
	def on_modified(self, event):
		if event.is_directory:
			return
		if event.src_path.endswith('main.md'):
			print os.system("pandoc -s main.md -o main.pdf")

if __name__ == "__main__":
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