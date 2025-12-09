Python 3.11.0 (main, Oct 24 2022, 18:26:48) [MSC v.1933 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license" for more information.
>>> from cryptography.fernet import Fernet
>>> key = Fernet.generate_key()
>>> key
b'ZB5Wbw010sEszaHDSNXFdiiwt5AHvJRqshD_exNSH8c='
>>> f = Fernet(key)
>>> token = f.encrypt(b"TEST!!")
>>> token
b'gAAAAABpOGpcc67-SaJBPXBDRl0sv7RXL9rvr6WRORUluqPJtVO0vUGBLD3ea-fXHdYlZcU3UcEj7pWf1gUqNWz0tzaO6rLpgA=='
>>> f.decrypt(token)
b'TEST!!'
>>>
