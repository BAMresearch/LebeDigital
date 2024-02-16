from server import app, db, User

# run this to add a user to db

username = ""
plain_password = ""

with app.app_context():
    # create new user
    new_user = User(username=username)
    new_user.set_password(plain_password)

    # Add to db
    db.session.add(new_user)
    db.session.commit()

print(f"Benutzer {username} wurde erfolgreich hinzugefügt.")
